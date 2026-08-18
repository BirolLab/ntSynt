#!/usr/bin/env python3
"""
Merge collinear synteny blocks from an ntSynt synteny block TSV file.

TSV format (tab-separated):
  block_id  assembly  contig  start  end  strand  num_minimizers  broken_reason

Merging criteria (blocks are NOT merged if any of the following are true):
  - Contig IDs differ between consecutive blocks (id_change)
  - Orientations differ between consecutive blocks (ori_change)
  - Any gap difference is negative, i.e. inconsistent order (inconsistent_order)
  - The gap difference across assemblies exceeds --indel tolerance (indel)
  - The maximum gap in any assembly meets or exceeds --merge threshold (merge)
"""

import argparse
import sys
from dataclasses import dataclass, field
from typing import Optional


@dataclass
class AssemblyBlock:
    "Represents the coordinates and orientation of a synteny block in a specific assembly."
    assembly: str
    contig_id: str
    start: int
    end: int
    strand: str


@dataclass
class SyntenyBlock:
    "Represents synteny block coordinates, and reason for discontinuity with previous synteny block"
    block_id: int
    assembly_blocks: dict = field(default_factory=dict)  # assembly -> AssemblyBlock
    broken_reason: Optional[str] = None


def get_difference_between_blocks(curr: AssemblyBlock, nxt: AssemblyBlock) -> int:
    """
    Compute the gap between two consecutive assembly blocks, respecting orientation.
    For '+' strand: gap = next.start - curr.end
    For '-' strand: gap = curr.start - next.end
    """
    if curr.strand == "+":
        return nxt.start - curr.end
    return curr.start - nxt.end


def parse_blocks(filepath: str) -> list[SyntenyBlock]:
    """Parse the TSV file into a list of SyntenyBlock objects, in order."""
    blocks: dict[int, SyntenyBlock] = {}
    first_assembly: Optional[str] = None

    with open(filepath, 'r', encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 6:
                raise ValueError(f"Expected at least 6 tab-separated columns, got {len(parts)}: {line!r}")

            block_id   = int(parts[0])
            assembly   = parts[1]
            contig_id  = parts[2]
            start      = int(parts[3])
            end        = int(parts[4])
            strand     = parts[5]

            if first_assembly is None:
                first_assembly = assembly

            if block_id not in blocks:
                blocks[block_id] = SyntenyBlock(block_id=block_id)

            blocks[block_id].assembly_blocks[assembly] = AssemblyBlock(
                assembly=assembly,
                contig_id=contig_id,
                start=start,
                end=end,
                strand=strand,
            )

    return [block for _, block in sorted(blocks.items(), key=lambda x: x[0])], first_assembly

def sort_blocks(blocks: list[SyntenyBlock], reference_assembly: str) -> list[SyntenyBlock]:
    """
    Sort blocks by the contig and start coordinate of the reference assembly
    (the first assembly encountered in the file).
    After sorting, the broken_reason of the first block is reset to None,
    since it no longer has a predecessor.
    """
    sorted_blocks = sorted(
        blocks,
        key=lambda b: (
            b.assembly_blocks[reference_assembly].contig_id,
            b.assembly_blocks[reference_assembly].start,
        )
    )
    # Reset broken_reason — will be fully recomputed during merging
    for block in sorted_blocks:
        block.broken_reason = None
    return sorted_blocks

def merge_collinear_blocks(blocks: list[SyntenyBlock],
                           indel_tolerance: int,
                           merge_threshold: int) -> list[SyntenyBlock]:
    """
    Merge consecutive collinear blocks.

    Two blocks are merged unless:
      - contig IDs differ                      -> broken_reason = "id_change"
      - orientations differ                    -> broken_reason = "ori_change"
      - unexpected order given orientations    -> broken_reason = "inconsistent_order"
      - max(gap) - min(gap) > indel_tolerance  -> broken_reason = "indel"
      - max(gap) >= merge_threshold            -> broken_reason = "merge"

    When merging, block coordinates are extended respecting strand:
      '+': curr.end   = next.end    (extend rightward)
      '-': curr.start = next.start  (extend leftward)
    """
    if not blocks:
        return []

    out_blocks: list[SyntenyBlock] = []
    curr_block = blocks[0]

    for block in blocks[1:]:
        orientations_match = True
        contig_id_match    = True
        differences        = []

        for assembly, curr_ab in curr_block.assembly_blocks.items():
            next_ab = block.assembly_blocks[assembly]

            if curr_ab.strand != next_ab.strand:
                orientations_match = False
            if curr_ab.contig_id != next_ab.contig_id:
                contig_id_match = False

            differences.append(get_difference_between_blocks(curr_ab, next_ab))

        gap_range = max(differences) - min(differences)
        max_gap   = max(differences)

        # Determine if we should break here, and why
        broken_reason = None
        if not contig_id_match:
            broken_reason = "id_change"
        elif not orientations_match:
            broken_reason = "ori_change"
        elif any(diff < 0 for diff in differences):
            broken_reason = "inconsistent_order"
        elif gap_range > indel_tolerance:
            broken_reason = "indel"
        elif max_gap >= merge_threshold:
            broken_reason = "merge"

        if broken_reason is not None:
            # Do not merge, start a new block
            block.broken_reason = broken_reason
            out_blocks.append(curr_block)
            curr_block = block
        else:
            # Merge: extend curr_block coordinates to absorb block
            for assembly, curr_ab in curr_block.assembly_blocks.items():
                next_ab = block.assembly_blocks[assembly]
                if curr_ab.strand == "+":
                    curr_ab.end = next_ab.end
                else:
                    curr_ab.start = next_ab.start

    out_blocks.append(curr_block)
    return out_blocks


def write_blocks(blocks: list[SyntenyBlock]) -> None:
    """Write merged blocks to TSV, renumbering block IDs from 0."""
    for new_id, block in enumerate(blocks):
        broken = block.broken_reason if block.broken_reason is not None else "None"
        for _, ab in block.assembly_blocks.items():
            print(str(new_id),
                ab.assembly,
                ab.contig_id,
                str(ab.start),
                str(ab.end),
                ab.strand,
                "0",
                broken, sep="\t", file=sys.stdout)


def main():
    "Parse arguments and run merging of collinear synteny blocks."
    parser = argparse.ArgumentParser(
        description="Merge collinear synteny blocks from an ntSynt TSV file."
    )
    parser.add_argument("--tsv", help="Input synteny blocks TSV file", required=True)
    parser.add_argument(
        "--indel",
        type=int,
        default=50000,
        help="Indel size threshold (bp) [%(default)s]"
    )
    parser.add_argument(
        "--merge",
        type=int,
        default=1000000,
        help="Maximum distance between collinear synteny blocks for merging (bp) [%(default)s]"
    )
    args = parser.parse_args()

    blocks, first_assembly = parse_blocks(args.tsv)
    blocks = sort_blocks(blocks, first_assembly)
    merged = merge_collinear_blocks(blocks, args.indel, args.merge)
    write_blocks(merged)


if __name__ == "__main__":
    main()
