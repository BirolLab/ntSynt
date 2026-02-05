#!/usr/bin/env python3
"""
ntSynt: Identifying synteny between genome assemblies using minimizer graphs
Written by Lauren Coombe @lcoombe
"""

from collections import defaultdict
import copy
import datetime
import os
import re
import shlex
import subprocess
import sys
import intervaltree
import pybedtools
import btllib
import ncls
import ntjoin_utils
import ntjoin
import igraph as ig
from synteny_block import SyntenyBlock


# Regexes
fai_re = re.compile(r'^(\S+).k\d+.w\d+.tsv')

class NtSyntSynteny(ntjoin.Ntjoin):
    "Instance for ntSynt synteny detection"

    def __init__(self, args):
        super().__init__(args)
        self.weights_list = [1] * len(self.args.FILES)
        self.args.t = 1
        self.args.FILES = sorted(self.args.FILES, reverse=True)
        self.print_parameters_synteny()
        self.repeat_bf = None
        if collinear_match := re.search(r"^(\d+)w$", self.args.collinear_merge):
            self.args.collinear_merge = int(collinear_match.group(1)) * self.args.w
        elif collinear_match := re.search(r"^(\d+)$", self.args.collinear_merge):
            self.args.collinear_merge = int(collinear_match.group(1))
        else:
            raise ValueError("--collinear-merge must be provided with an integer value or string in the form '<num>w'")

    def print_parameters_synteny(self):
        "Pring the set parameters for the ntSynt synteny run"
        if self.args.n == 0:
            self.args.n = len(self.args.FILES)
        print("Parameters:")
        print("\tMinimizer TSV files: ", self.args.FILES)
        print("\t--fastas", self.args.fastas)
        print("\t-n", self.args.n)
        print("\t-p", self.args.p)
        print("\t-k", self.args.k)
        print("\t-w", self.args.w)
        print("\t--btllib_t", self.args.btllib_t)
        print("\t--w-rounds", self.args.w_rounds)
        print("\t-m", self.args.m)
        print("\t-z", self.args.z)
        print("\t--collinear-merge", self.args.collinear_merge, flush=True)
        if self.args.common:
            print("\t--common", self.args.common, flush=True)
        if self.args.repeat:
            print("\t--repeat", self.args.repeat, flush=True)


    def find_synteny_blocks(self, path):
        "Given a path (sequence of mx), print the order/orientation/regions of contigs for an assembly"
        out_blocks = []  # List of SyntenyBlock
        to_remove_nodes = [] # List of nodes to remove - if block isn't oriented
        prelim_blocks = SyntenyBlock(self.args.k, self.args.m, *list(self.list_mx_info.keys()))
        past_start_flag = False
        for mx in path:
            if prelim_blocks.continue_block(mx, self.list_mx_info):
                prelim_blocks.extend_block(mx, self.list_mx_info)
            else:
                # This is either the first mx, or we are past a stretch of repeating contigs
                if past_start_flag:
                    prelim_blocks.determine_orientations()
                    if prelim_blocks.all_oriented():
                        out_blocks.append(prelim_blocks)
                    else:
                        if self.args.dev:
                            print("Not oriented: ", prelim_blocks, flush=True)
                        to_remove_nodes.extend([ntjoin_utils.vertex_index(self.graph, mx.mx)
                                        for mx in prelim_blocks.assembly_blocks[
                                            list(prelim_blocks.assembly_blocks.keys()).pop()].minimizers])
                prelim_blocks = SyntenyBlock(self.args.k, self.args.m, *list(self.list_mx_info.keys()))
                prelim_blocks.start_block(mx, self.list_mx_info)

        prelim_blocks.determine_orientations()
        if prelim_blocks.all_oriented():
            out_blocks.append(prelim_blocks)
        else:
            if self.args.dev:
                print("Not oriented: ", prelim_blocks, flush=True)
            to_remove_nodes.extend([ntjoin_utils.vertex_index(self.graph, mx.mx)
                                    for mx in prelim_blocks.assembly_blocks[
                                        list(prelim_blocks.assembly_blocks.keys()).pop()].minimizers])

        # Remove the nodes from unoriented blocks from the graph
        if to_remove_nodes:
            new_graph = self.graph.copy()
            new_graph.delete_vertices(to_remove_nodes)
            self.graph = new_graph

        return out_blocks

    @staticmethod
    def find_fa_name(assembly_mx_name):
        "Given the mx file name, return the corresponding fai file name"
        if fai_match := re.search(fai_re, assembly_mx_name):
            return f"{fai_match.group(1)}"
        print("ERROR: Target assembly minimizer TSV file must follow the naming convention:")
        print("\ttarget_assembly.fa.k<k>.w<w>.tsv, where <k> and <w> are parameters used for minimizering")
        sys.exit(1)

    @staticmethod
    def get_synteny_bed_lists(paths):
        "Given a set of synteny blocks, return a dictionary with a BED interval list per contig, per assembly"
        synteny_beds = {}
        for block in paths:
            for assembly, assembly_block in block.assembly_blocks.items():
                if assembly not in synteny_beds:
                    synteny_beds[assembly] = {}
                if assembly_block.contig_id not in synteny_beds[assembly]:
                    synteny_beds[assembly][assembly_block.contig_id] = []
                synteny_beds[assembly][assembly_block.contig_id].append(
                    ntjoin_utils.Bed(assembly_block.contig_id,
                                    assembly_block.get_block_start(),
                                    assembly_block.get_block_end()))

        return synteny_beds

    def mask_assemblies_with_synteny_extents(self, synteny_beds, w):
        "Mask each reference assembly with determined synteny blocks"
        mx_to_fa_dict = {}
        assembly_to_fastas = {os.path.basename(fasta): fasta for fasta in self.args.fastas}
        for assembly, contig_dict in synteny_beds.items():
            bed_str = [f"{ctg}\t{bed.start}\t{bed.end}\tSYNTENY" for ctg in contig_dict \
                        for bed in contig_dict[ctg] if bed.end - bed.start > max(2*w, w+self.args.k+1)]
            bed_str = "\n".join(bed_str)
            fa_filename = self.find_fa_name(assembly)
            fa_filename_full = assembly_to_fastas[fa_filename]
            synteny_bed = pybedtools.BedTool(bed_str, from_string=True).slop(g=f"{fa_filename}.fai",
                                                                             l=-1*(w+self.args.k),
                                                                             r=-1*(w+self.args.k)).sort()
            synteny_bed.mask_fasta(fi=fa_filename_full, fo=f"{fa_filename}_masked.fa.tmp")

            # pybedtools may output multi-line fasta which breaks btllib reading, need seqtk to make single line
            with open(f"{fa_filename}_masked.fa", "w", encoding="utf-8") as fout:
                process = subprocess.run(shlex.split(f"seqtk seq {fa_filename}_masked.fa.tmp"),
                                        stdout=fout,
                                        check=True, text=True)
                assert process.returncode == 0

            mx_to_fa_dict[assembly] = f"{fa_filename}_masked.fa"
        return mx_to_fa_dict

    @staticmethod
    def delete_w_iteration_files(*filenames):
        "Delete the given files for the specific lower w iteration"
        for filename in filenames:
            cmd = shlex.split(f"rm {filename}")
            ret_code = subprocess.call(cmd)
            assert ret_code == 0

    def generate_new_minimizers(self, tsv_to_fa_dict, w, retain_files=False):
        "Given the masked fasta files, generate minimizers at new w for each"
        list_mxs = {}
        new_list_mxs_info = {}
        for assembly_tsv, assembly_masked in tsv_to_fa_dict.items():
            if self.args.filter == "Indexlr" and not self.args.common:
                indexlr_filename = ntjoin_utils.run_indexlr(assembly_masked, self.args.k, w,
                                                             self.args.btllib_t, r=self.args.repeat)
            elif self.args.common and self.args.filter != "Indexlr":
                indexlr_filename = ntjoin_utils.run_indexlr(assembly_masked, self.args.k, w,
                                                             self.args.btllib_t, s=self.args.common)
            elif self.args.common and self.args.filter == "Indexlr":
                indexlr_filename = ntjoin_utils.run_indexlr(assembly_masked, self.args.k, w,
                                                             self.args.btllib_t, s=self.args.common, r=self.args.repeat)
            else:
                indexlr_filename = ntjoin_utils.run_indexlr(assembly_masked, self.args.k, w, self.args.btllib_t)

            if self.args.filter == "Filter":
                mx_info, mxs_filt = ntjoin_utils.read_minimizers(indexlr_filename, self.repeat_bf)
            else:
                mx_info, mxs_filt = ntjoin_utils.read_minimizers(indexlr_filename)
            new_list_mxs_info[assembly_tsv] = mx_info
            list_mxs[assembly_tsv] = mxs_filt
            if not retain_files:
                self.delete_w_iteration_files(indexlr_filename, assembly_masked, f"{assembly_masked}.tmp")
        return list_mxs, new_list_mxs_info

    @staticmethod
    def update_intervals(assembly_name, ctg, mx1, mx2, intervals):
        "Update the given dictionary of trees with the new extent"
        start_pos = min(mx1.position, mx2.position)
        end_pos = max(mx1.position, mx2.position)
        if end_pos - start_pos < 2: # If the interval is too short, will generate an error if try to add the interval
            return
        if assembly_name not in intervals or ctg not in intervals[assembly_name]:
            intervals[assembly_name][ctg] = []
        intervals[assembly_name][ctg].append((start_pos+1, end_pos, 1))

    def find_mx_in_blocks(self, paths):
        "Given the synteny blocks, find the minimizers at the terminal ends of each block, and internal mxs"
        terminal_mxs = set()
        internal_mxs = set()
        intervals = defaultdict(dict)

        for block in paths:
            curr_mx_len = len(terminal_mxs)
            for assembly, assembly_block in block.assembly_blocks.items():
                contig, mx1, mx2 = assembly_block.get_block_terminal_mx()
                terminal_mxs.add(mx1.mx)
                terminal_mxs.add(mx2.mx)
                self.update_intervals(assembly, contig, mx1, mx2, intervals)
                internal = assembly_block.get_block_internal_mx_hashes()
                internal_mxs.update(internal)
            assert len(terminal_mxs) == (curr_mx_len + 2)

        for assembly in intervals:
            for ctg in intervals[assembly]:
                start, end, data = zip(*intervals[assembly][ctg])
                intervals[assembly][ctg] = ncls.NCLS(list(start), list(end), list(data))
        return terminal_mxs, internal_mxs, intervals

    def get_overlapping_region(self, start: int, end: int, interval: intervaltree.Interval) -> intervaltree.Interval:
        "Given a start/end and list of intervals, return a list of intervals describing the overlapping region"
        new_start = max(start, interval.begin)
        new_end = min(end, interval.end)
        return intervaltree.Interval(new_start, new_end)

    def check_non_overlapping(self, paths):
        "Given the paths, do final check to ensure intervals are not overlapping, will print warnings if so"
        intervaltrees = defaultdict(dict) # assembly -> contig -> IntervalTree of synteny block extents
        for block in paths:
            for assembly, assembly_block in block.assembly_blocks.items():
                contig, start, end = assembly_block.get_block_contig_start_end()
                if assembly not in intervaltrees or contig not in intervaltrees[assembly]:
                    intervaltrees[assembly][contig] = intervaltree.IntervalTree()
                if not (all(asm_block.get_block_length() >= self.args.z # Don't worry about short ones
                               for _, asm_block in block.assembly_blocks.items())):
                    continue
                if hit_intervals := intervaltrees[assembly][contig][start:end]: # Check doesn't overlap with anything
                    for hit_interval in hit_intervals:
                        overlapping_region = self.get_overlapping_region(start, end, hit_interval)
                        if overlapping_region.end - overlapping_region.begin >= self.args.z:
                            print("WARNING: detected overlapping segments for this block:", assembly,
                                  contig, start, end,
                                    "\n", file=sys.stderr, flush=True)
                            break
                intervaltrees[assembly][contig][start:end] = (start, end)


    @staticmethod
    def filter_minimizers_synteny_blocks(list_mxs, black_list, list_mx_info, intervals):
        "Filter minimizers found in the mx black list"
        return_mxs = {}
        for assembly in list_mxs:
            assembly_mxs_filtered = []
            for mx_list in list_mxs[assembly]:
                new_list = []
                for mx in mx_list:
                    ctg, pos = list_mx_info[assembly][mx]
                    asm_ctg_intervals = intervals[assembly][ctg] if ctg in intervals[assembly] else None
                    if new_list and asm_ctg_intervals is not None:
                        prev_pos = list_mx_info[assembly][new_list[-1]][1]
                        start = min(prev_pos, pos)
                        end = max(prev_pos, pos)
                        if asm_ctg_intervals.has_overlap(start, end):
                            assembly_mxs_filtered.append(new_list)
                            new_list = []
                    if mx not in black_list and (asm_ctg_intervals is None or
                                                 not asm_ctg_intervals.has_overlap(pos, pos+1)):
                        new_list.append(mx)
                assembly_mxs_filtered.append(new_list)

            return_mxs[assembly] = assembly_mxs_filtered
        return return_mxs

    @staticmethod
    def update_list_mx_info(list_mxs, list_mx_info, new_list_mx_info):
        "Update the directory containing mx -> contig, position associations"
        valid_mxs = set({mx for _, list_mx_val in list_mxs.items() \
                        for list_mx in list_mx_val for mx in list_mx})
        for assembly, mx_dict in new_list_mx_info.items():
            for mx in mx_dict:
                if mx in valid_mxs: #and mx not in list_mx_info[assembly]:
                    list_mx_info[assembly][mx] = mx_dict[mx]

    def filter_graph_global_flag_overlaps(self, graph):
        "Filter the graph globally based on minimum edge weight, flagging incident nodes"
        print(datetime.datetime.today(), ": Filtering the graph", file=sys.stdout, flush=True)
        flagged_node_pairs = []
        to_remove_edges = []
        for edge in graph.es():
            if edge['weight'] < self.args.n:
                to_remove_edges.append(edge.index)
                flagged_node_pairs.append((edge.source, edge.target))
        new_graph = graph.copy()
        new_graph.delete_edges(to_remove_edges)
        return new_graph, flagged_node_pairs

    def has_overlap(self, source, target):
        "Return True if the positions supplied are less than k apart"
        for _, mx_info in self.list_mx_info.items():
            if abs(mx_info[source][1] - mx_info[target][1]) < self.args.k:
                return True
        return False

    def erode_edges(self, source, target):
        "Given a source and target, erode incident edges as needed to ensure that the vertices do not overlap"
        erode_target = True
        curr_source, curr_target = source, target
        return_edges, visited = set(), set([curr_source, curr_target])
        source_name = ntjoin_utils.vertex_name(self.graph, source)
        target_name = ntjoin_utils.vertex_name(self.graph, target)

        while self.has_overlap(source_name, target_name):
            erode_vertex = curr_target if erode_target else curr_source
            for edge in self.graph.vs()[erode_vertex].incident():
                return_edges.add(edge)
            possible_neighbours = [vertex for vertex in self.graph.vs()[erode_vertex].neighbors() \
                                    if vertex.index not in visited]
            if len(possible_neighbours) == 0:
                break
            assert len(possible_neighbours) == 1
            if erode_target:
                curr_target = possible_neighbours.pop().index
                target_name = ntjoin_utils.vertex_name(self.graph, curr_target)
                erode_target = False
                visited.add(curr_target)
            else:
                curr_source = possible_neighbours.pop().index
                source_name = ntjoin_utils.vertex_name(self.graph, curr_source)
                erode_target = True
                visited.add(curr_source)

        return list(return_edges)


    def refine_graph(self, flagged_node_pairs):
        "Refine the undirected minimizer graph by filtering edges between blocks that overlap"
        to_remove_edges = []
        if not flagged_node_pairs:
            return self.graph

        for source, target in flagged_node_pairs:
            # normalize source/target to ensure deterministic
            if ntjoin_utils.vertex_name(self.graph, source) > ntjoin_utils.vertex_name(self.graph, target):
                source, target = target, source
            # Only consider starting the trimming if these are terminal nodes
            if self.graph.vs()[source].degree() != 1 or self.graph.vs()[target].degree() != 1:
                continue
            to_remove_edges.extend(self.erode_edges(source, target))

        new_graph = self.graph.copy()
        if not to_remove_edges:
            return self.graph
        new_graph.delete_edges(to_remove_edges)
        return new_graph

    @staticmethod
    def max_difference(mx_list1, mx_list2):
        "Return the maximum difference in interarrival distances betweens the minimizers in the list"
        interarrivals = [abs(pos1 - mx_list2.positions[i]) for i, pos1 in enumerate(mx_list1.positions)]
        return max(interarrivals) - min(interarrivals)

    def break_synteny_block(self, block, break_positions):
        "Given a synteny block, break at the specified minimizer positions"
        if not break_positions:
            return [block]
        return_blocks = []
        break_intervals = intervaltree.IntervalTree()

        break_intervals[0:block.get_number_of_minimizers()] = 0
        for pos in break_positions:
            break_intervals.slice(pos)
        for new_block_extents in sorted(break_intervals):
            new_synteny_block = SyntenyBlock(self.args.k, self.args.m, *list(block.assembly_blocks.keys()))
            for assembly in block.assembly_blocks:
                new_assembly_block = copy.copy(block.assembly_blocks[assembly])
                new_assembly_block.minimizers = \
                    block.assembly_blocks[assembly].minimizers[new_block_extents.begin: new_block_extents.end]
                new_synteny_block.assign_block(assembly, new_assembly_block)
            return_blocks.append(new_synteny_block)
        return return_blocks


    def check_for_indels(self, paths):
        "Given a set of paths, check each for any detected indels, breaking paths if needed"
        return_paths = []
        remove_edges = []
        for block in paths:
            break_positions = []
            for i in range(block.get_number_of_minimizers() - 1):
                minimizers_to_compare = (block.get_node(i), block.get_node(i+1))
                if self.max_difference(*minimizers_to_compare) > self.args.bp:
                    break_positions.append(i+1)
                    mx_i, mx_i_next = minimizers_to_compare[0].mx, minimizers_to_compare[1].mx
                    remove_edges.append(ntjoin_utils.edge_index(self.graph, mx_i, mx_i_next))
            if not break_positions:
                return_paths.append(block)
            else:
                return_paths.extend(self.break_synteny_block(block, break_positions))
        self.graph = ntjoin_utils.remove_flagged_edges(self.graph, remove_edges)

        return return_paths

    def filter_synteny_blocks(self, paths, mx_threshold=1):
        "Filter out the synteny blocks that are comprised of less than mx_threshold minimizers"
        return_blocks = []
        to_remove_nodes = []
        for block in paths:
            if all(len(asm_block.minimizers) >= mx_threshold
                   for _, asm_block in block.assembly_blocks.items()):
                return_blocks.append(block)
            else:
                to_remove_nodes.extend([ntjoin_utils.vertex_index(self.graph, mx.mx)
                                        for mx in block.assembly_blocks[
                                            list(block.assembly_blocks.keys()).pop()].minimizers])
        new_graph = self.graph.copy()
        new_graph.delete_vertices(to_remove_nodes)
        self.graph = new_graph
        return return_blocks

    def get_difference_between_blocks(self, block1, block2):
        "Returns the gap between the blocks on the assembly"
        if block1.ori == "-" and block2.ori == "-":
            return  block1.get_block_start() - block2.get_block_end()
        return block2.get_block_start() - block1.get_block_end()

    def merge_collinear_blocks(self, blocks):
        "Merge collinear blocks that are less than threshold apart, and are not indels"
        out_blocks = []
        curr_block = blocks[0]
        for block in blocks[1:]:
            orientations = True # Are the orientations all the same?
            contig_id = True # Are the contig IDs all the same?
            differences = []
            for assembly, assembly_block in curr_block.assembly_blocks.items():
                # Tally orientation info
                orientations = False if assembly_block.ori != block.assembly_blocks[assembly].ori \
                    else orientations
                contig_id = False if assembly_block.contig_id != block.assembly_blocks[assembly].contig_id \
                    else contig_id
                differences.append(self.get_difference_between_blocks(assembly_block,
                                                                      block.assembly_blocks[assembly]))
            if not orientations or not contig_id or \
                (max(differences) - min(differences) > self.args.bp - self.args.k) or \
                    max(differences) >= self.args.collinear_merge:
                if not contig_id:
                    block.broken_reason = "id_change"
                elif not orientations:
                    block.broken_reason = "ori_change"
                elif any(diff < 0 for diff in differences):
                    block.broken_reason = "inconsistent_order"
                elif max(differences) - min(differences) > self.args.bp - self.args.k:
                    block.broken_reason = "indel"
                elif max(differences) >= self.args.collinear_merge:
                    block.broken_reason = "merge"

                out_blocks.append(curr_block)
                curr_block = block
            else:
                # Extend this block
                for assembly, assembly_block in block.assembly_blocks.items():
                    curr_block.assembly_blocks[assembly].minimizers.extend(assembly_block.minimizers)

        out_blocks.append(curr_block)
        return out_blocks



    def refine_block_coordinates(self, paths):
        "Ready to start refining the synteny block coordinates"
        prev_w = self.args.w
        for new_w in self.args.w_rounds:
            print(datetime.datetime.today(), ": Extending synteny blocks with w =", new_w, file=sys.stdout, flush=True)
            new_list_mxs, terminal_mxs = self.generate_additional_minimizers(
                    paths, new_w, prev_w, self.list_mx_info, self.args.dev)
            graph = ntjoin_utils.build_graph(new_list_mxs, self.weights, graph=self.graph, black_list=terminal_mxs)
            if self.args.simplify_graph:
                self.graph = self.run_graph_simplification(self.graph)
            # Do further graph refinement on last round only
            if new_w == self.args.w_rounds[-1]:
                self.graph, node_pairs = self.filter_graph_global_flag_overlaps(graph)
                self.graph = self.refine_graph(node_pairs)
            else:
                self.graph = self.filter_graph_global(graph)
            paths = self.find_paths_synteny_blocks(self.find_paths())
            paths = self.check_for_indels(paths)
            paths = self.filter_synteny_blocks(paths, 4) # TODO: magic number
            paths_sorted_for_printing = sorted(paths)
            with open(f"{self.args.p}.pre-collinear-merge.synteny_blocks.tsv", 'w', encoding="utf-8") as outfile:
                block_num = 0
                for block in paths_sorted_for_printing:
                    if not all(asm_block.get_block_length() >= self.args.z
                               for _, asm_block in block.assembly_blocks.items()):
                        continue
                    outfile.write(block.get_block_string(block_num))
                    block_num += 1
            if new_w == self.args.w_rounds[-1]:
                paths_sorted_for_printing_collinear = self.merge_collinear_blocks(paths_sorted_for_printing)
                # Filter by length, and do another round of merging
                paths_sorted_for_printing_collinear = [block for block in paths_sorted_for_printing_collinear \
                                                        if all(asm_block.get_block_length() >= self.args.z \
                                                            for _, asm_block in block.assembly_blocks.items())]
                paths_sorted_for_printing_collinear = self.merge_collinear_blocks(paths_sorted_for_printing_collinear)

                # Check for non-overlapping in the last round if --dev
                if self.args.dev:
                    self.check_non_overlapping(paths_sorted_for_printing_collinear)

                with open(f"{self.args.p}.synteny_blocks.tsv", 'w', encoding="utf-8") as outfile:
                    block_num = 0
                    for block in paths_sorted_for_printing_collinear:
                        if not all(asm_block.get_block_length() >= self.args.z
                                for _, asm_block in block.assembly_blocks.items()):
                            continue
                        outfile.write(block.get_block_string(block_num, verbose=True))
                        block_num += 1

            prev_w = new_w


        print(datetime.datetime.today(), ": Done extended synteny blocks", file=sys.stdout, flush=True)
        print(datetime.datetime.today(),
                f": Final synteny blocks can be found in: {self.args.p}.synteny_blocks.tsv", flush=True)

    def generate_additional_minimizers(self, paths, new_w, prev_w, list_mx_info, dev=False):
        "Given the existing synteny blocks, generate minimizers for increased block resolution"
        synteny_beds = self.get_synteny_bed_lists(paths)
        mx_to_fa_dict = self.mask_assemblies_with_synteny_extents(synteny_beds, prev_w)
        list_mxs, new_list_mx_info = self.generate_new_minimizers(mx_to_fa_dict, new_w, retain_files=dev)
        terminal_mx, internal_mx, intervals = self.find_mx_in_blocks(paths)
        list_mxs = self.filter_minimizers_synteny_blocks(list_mxs, internal_mx, new_list_mx_info, intervals)
        list_mxs = ntjoin_utils.filter_minimizers(list_mxs) # Filter for mx in all assemblies
        self.update_list_mx_info(list_mxs, list_mx_info, new_list_mx_info)
        return list_mxs, terminal_mx

    def find_paths_synteny_blocks(self, paths):
        "Given a list of paths, return a list of representative synteny blocks"
        print(datetime.datetime.today(), ": Finding synteny blocks", flush=True)
        return [block for path in paths for blocks, _ in path for block in self.find_synteny_blocks(blocks)]

    @staticmethod
    def node_partially_anchored(graph, vertex_id, max_edge_weight):
        '''Return true if the given node is partially anchored. Partially anchored is defined as a node
        that has only one incident edge with max (num_assemblies) weight'''
        incident_weights = [e["weight"] for e in graph.es()[graph.incident(vertex_id)]]
        if incident_weights.count(max_edge_weight) == 1:
            return True
        return False

    def print_interarrivals(self, paths):
        "Print the interarrival distances between nodes"
        with open(f"{self.args.p}.interarrivals.tsv", 'w', encoding="utf-8") as fout:
            for block in paths:
                for _, assembly_block in block.assembly_blocks.items():
                    for mx1, mx2 in zip(assembly_block.minimizers,
                                        assembly_block.minimizers[1:]):
                        fout.write(f"{abs(mx2.position - mx1.position)}\n")

    def run_graph_simplification(self, graph):
        '''Run graph simplification.
        For each edge between partially anchored nodes, check if there is another simple, 2 step path
        If so, delete the middle node (likely problematic), and set the transitive edge to the max weight'''
        print(datetime.datetime.today(), ": Running graph simplificaton", file=sys.stdout, flush=True)
        max_edge_weight = sum(self.weights.values())
        to_remove_nodes = []

        # get_all_simple_paths API changed in igraph 1.0.0, account for that here
        max_path_len = {"cutoff": 2} if ig.__version__ < "1.0.0" else {"maxlen": 2}

        for edge in graph.es():
            source, target = edge.source, edge.target
            # They are potentially simple bubble-like structures, with anchors on either side,
            # since total degree for all nodes must be 2*number_input_assemblies
            if graph.vs()[source].degree() == 3 and \
                graph.vs()[target].degree() == 3 and \
                self.node_partially_anchored(graph, source, max_edge_weight) and \
                self.node_partially_anchored(graph, target, max_edge_weight):
                paths = graph.get_all_simple_paths(source, target, **max_path_len)
                if len(paths) == 2: # Found the one with this edge plus one more which is 3 nodes long
                    for path in paths:
                        if len(path) == 3: # This is the one with the node we want to get rid of
                            to_remove_nodes.append(path[1]) # Appending the middle node
                            edge["weight"] = max_edge_weight

        new_graph = graph.copy()
        new_graph.delete_vertices(to_remove_nodes)
        return new_graph


    def main_synteny(self):
        "Run the steps for ntSynt synteny detection"

        # Check that the w_rounds specified don't have any duplicates
        if len(self.args.w_rounds) != len(set(self.args.w_rounds)):
            print("Error: duplicate values found in w_rounds!", file=sys.stderr, flush=True)
            sys.exit(1)

        if self.args.filter and not self.args.repeat:
            raise ValueError("If --filter is specified, must supply repeat Bloom filter with --repeat")

        # Load the minimizers
        if self.args.filter == "Filter":
            self.repeat_bf = btllib.KmerBloomFilter(self.args.repeat)
            self.load_minimizers(self.repeat_bf)
        else:
            self.load_minimizers()

        # Make the minimizer graph
        self.make_minimizer_graph()

        # Simplify the minimizer graph, then filter
        if self.args.simplify_graph:
            self.graph = self.run_graph_simplification(self.graph)
        self.graph = self.filter_graph_global(self.graph)

        # Find the initial paths
        paths = self.ntjoin_find_paths()

        # Detect the synteny blocks from graph paths, refine them
        paths = self.find_paths_synteny_blocks(paths)
        paths = self.check_for_indels(paths)
        paths = self.filter_synteny_blocks(paths, 4)
        if self.args.interarrivals:
            self.print_interarrivals(paths)
        paths_sorted_for_printing = sorted(paths)

        if not paths_sorted_for_printing:
            print("Error - no paths found. Try adjusting the specified k/w parameters.")
            sys.exit(1)

        with open(f"{self.args.p}.synteny_blocks.tsv", 'w', encoding="utf-8") as outfile:
            block_num = 0
            for block in paths_sorted_for_printing:
                if not all(asm_block.get_block_length() >= self.args.z
                           for _, asm_block in block.assembly_blocks.items()):
                    continue
                outfile.write(block.get_block_string(block_num))
                block_num += 1
        print(datetime.datetime.today(), ": Done initial synteny blocks", file=sys.stdout, flush=True)

        # Further refine the synteny blocks
        self.refine_block_coordinates(paths)

        print(datetime.datetime.today(), ": DONE!", file=sys.stdout, flush=True)
