#include "btllib/bloom_filter.hpp"
#include "btllib/seq_reader.hpp"
#include "btllib/util.hpp"
#include <argparse/argparse.hpp>
#include <cmath>
#include <iostream>
#include <vector>

#if _OPENMP
#include <omp.h>
#endif

/*
Creating a common Bloom filter (with k-mers found in all input assemblies) using
a cascading Bloom filter approach
*/

// Converting bits to bytes
const unsigned NUM_BITS_PER_BYTE = 8;

/*
Approximate the required BF size based on the input genome size.
Using the formula found in Broder & Mitzenmacher, 2004, simplified for
one hash function.
*/
long long
approximate_bf_size(std::string genome_file, double fpr, int hashes, int threads)
{
  long long genome_size = 0;
  btllib::SeqReader reader(
    genome_file, btllib::SeqReader::Flag::LONG_MODE, threads);
  for (const auto record : reader) {
    genome_size += record.seq.length();
  }
  std::cout << "Genome size (bp): " << genome_size << std::endl;
  double m_bits = -1.0 * hashes * genome_size / log(1.0 - pow(fpr, 1.0 / hashes));

  long long size_bits = (long long) std::ceil(m_bits);
  return size_bits / NUM_BITS_PER_BYTE;
}

int
main(int argc, const char** argv)
{

  argparse::ArgumentParser parser("ntsynt_make_common_bf");
  parser.add_argument("--genome")
    .nargs(argparse::nargs_pattern::at_least_one)
    .help("Input genome file(s)")
    .required();
  parser.add_argument("-k")
    .help("k-mer size (bp)")
    .required()
    .scan<'u', unsigned>();

  parser.add_argument("--fpr")
    .help("False positive rate for Bloom filter")
    .default_value((double)0.025)
    .scan<'g', double>();

  parser.add_argument("-p")
    .help("Prefix for output Bloom filter")
    .default_value("common_bf");

  parser.add_argument("--bf")
    .help("Bloom filter size in bytes (optional)")
    .scan<'d', long long>();
  
  parser.add_argument("--hashes")
    .help("Number of hash functions")
    .default_value(1U)
    .scan<'u', unsigned>();

  parser.add_argument("-t")
    .help("Number of threads")
    .default_value(12U)
    .scan<'u', unsigned>();

  /* Parse the command-line arguments */
  try {
    parser.parse_args(argc, argv);
  } catch (const std::runtime_error& err) {
    std::cerr << err.what() << std::endl;
    std::cerr << parser;
    std::exit(1);
  }

  std::vector<std::string> genome_files =
    parser.get<std::vector<std::string>>("genome");
  unsigned num_threads = parser.get<unsigned>("t");
  double fpr = parser.get<double>("fpr");
  unsigned k = parser.get<unsigned>("k");
  unsigned hashes = parser.get<unsigned>("hashes");
  std::string prefix = parser.get<std::string>("p");

  std::cout << "Parameters:" << std::endl;
  std::cout << "\t\t--genome ";
  for (const auto& genome : genome_files) {
    std::cout << genome << " ";
  }
  std::cout << std::endl;
  std::cout << "\t\t-t " << num_threads << std::endl;
  std::cout << "\t\t-k " << k << std::endl;
  std::cout << "\t\t--fpr " << fpr << std::endl;
  std::cout << "\t\t--hashes " << hashes << std::endl;
  std::cout << "\t\t-p " << prefix << std::endl;

  #if _OPENMP
  omp_set_num_threads(num_threads);
  #endif

  /* Sort the vector of genomes, so that the output BF will be the same 
  even if files are specified in a different order */
  std::sort(genome_files.begin(), genome_files.end());

  /* Calculate BF size based on genome size of first specified genome */
  long long bf_size;
  if (parser.is_used("--bf")) {
    bf_size = parser.get<long long>("bf");
    std::cout << "\t\t--bf " << bf_size << std::endl;
  } else {
    std::cout << "Calculating BF size based on input genome size" << std::endl;
    bf_size = approximate_bf_size(genome_files[0], fpr, hashes, num_threads);
  }
  std::cout << "BF size (bytes): " << bf_size << std::endl;


  /* Load the initial level 1 BF with k-mers from the first genome  */
  btllib::KmerBloomFilter* bf =
    new btllib::KmerBloomFilter(bf_size, hashes, k);
  btllib::log_info("Reading " + genome_files[0]);
  btllib::SeqReader reader(
    genome_files[0], btllib::SeqReader::Flag::LONG_MODE, num_threads);

  #pragma omp parallel
  for (const auto record : reader) {
    bf->insert(record.seq);
  }
  std::cout << "Bloom filter FPR: " << bf->get_fpr() << std::endl;

  /* For each subsequent genome, only insert into next level BF if the k-mer is
  in the previous level BF. */
  btllib::KmerBloomFilter* new_bf =
    new btllib::KmerBloomFilter(bf_size, hashes, k);

  size_t num_assemblies = genome_files.size();
  for (size_t i = 1; i < num_assemblies; ++i) {
    std::string genome = genome_files[i];
    btllib::log_info("Reading " + genome);
    btllib::SeqReader reader(
      genome, btllib::SeqReader::Flag::LONG_MODE, num_threads);
    #pragma omp parallel
    for (const auto record : reader) {
      btllib::NtHash nthash(record.seq, hashes, k);
      while (nthash.roll()) {
        if (bf->contains(nthash.hashes())) {
          new_bf->insert(nthash.hashes());
        }
      }
    }
    std::cout << "Bloom filter FPR: " << new_bf->get_fpr() << std::endl;
    delete bf;
    bf = new_bf;
    if (i < num_assemblies - 1) {
      new_bf = new btllib::KmerBloomFilter(bf_size, hashes, k);
    }
  }

  std::cout << "Final Bloom filter FPR: " << bf->get_fpr() << std::endl;
  btllib::log_info("Saving Bloom filter");
  bf->save(prefix + ".bf");
  btllib::log_info("Done!");
  delete bf;
}