#include <seqan3/std/charconv>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/range/views/get.hpp>
#include <seqan3/range/views/persist.hpp>
#include <seqan3/range/views/move.hpp>
#include <seqan3/search/dream_index/technical_binning_directory.hpp>

struct cmd_arguments
{
    std::filesystem::path bin_path{};
    std::filesystem::path out_path{"./ibf.out"};
    uint64_t w{23};
    uint8_t k{20};
    uint64_t bins{64};
    uint64_t bits{4096};
    std::string size{};
    uint64_t hash{2};
    uint8_t threads{1};
    bool gz{false};
    bool bz2{false};
};

struct my_traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

void run_program(cmd_arguments & args)
{
    std::string extension{".fasta"};
    if (args.gz)
        extension += ".gz";
    if (args.bz2)
        extension += ".bz2";

    using sequence_file_t = seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq>>;

    auto technical_bins = std::views::iota(0u, args.bins) |
                          std::views::transform([&] (size_t const i) {
                                return sequence_file_t{args.bin_path / ("bin_" + std::to_string(i) + extension)} |
                                       seqan3::views::persist |
                                       seqan3::views::get<seqan3::field::seq> |
                                       seqan3::views::move;
                          });

    seqan3::ibf_config cfg{seqan3::bin_count{args.bins},
                           seqan3::bin_size{args.bits},
                           seqan3::hash_function_count{args.hash},
                           args.threads};

    seqan3::technical_binning_directory tbd{technical_bins, seqan3::views::kmer_hash(seqan3::ungapped{args.k}), cfg};

    std::ofstream os{args.out_path, std::ios::binary};
    cereal::BinaryOutputArchive oarchive{os};
    oarchive(tbd);
}

void initialize_argument_parser(seqan3::argument_parser & parser, cmd_arguments & args)
{
    parser.info.author = "Enrico Seiler";
    parser.info.author = "enrico.seiler@fu-berlin.de";
    parser.info.short_description = "Build an Interleaved Bloom Filter.";
    parser.info.version = "0.0.1";
    parser.add_option(args.bin_path, '\0', "input", "Please provide a path to a directory containing one FASTA file for"
                                                    " each bin.");
    parser.add_option(args.out_path, '\0', "output", "Please provide a valid output path.",
                      seqan3::option_spec::DEFAULT);
    parser.add_option(args.k, '\0', "kmer", "Choose the kmer size.", seqan3::option_spec::DEFAULT,
                      seqan3::arithmetic_range_validator{1, 32});
    parser.add_option(args.threads, '\0', "threads", "Choose the number of threads.", seqan3::option_spec::DEFAULT,
                      seqan3::arithmetic_range_validator{1, 2048});
    parser.add_option(args.bins, '\0', "bins", "Choose the number of bins.", seqan3::option_spec::REQUIRED,
                      seqan3::arithmetic_range_validator{1, 65536});
    parser.add_option(args.bits, '\0', "bits", "Choose the size in bits of one bin. Mutually exclusive with --size.",
                      seqan3::option_spec::DEFAULT, seqan3::arithmetic_range_validator{1, 35184372088832});
    parser.add_option(args.size, '\0', "size", "Choose the size of the resulting IBF. Mutually exclusive with --bits.");
    parser.add_option(args.hash, '\0', "hash", "Choose the number of hashes.", seqan3::option_spec::DEFAULT,
                      seqan3::arithmetic_range_validator{1, 4});
    parser.add_flag(args.gz, '\0', "gz", "Expect FASTA files to be gz compressed.");
    parser.add_flag(args.bz2, '\0', "bz2", "Expect FASTA files to be bz2 compressed.");
}

int main(int argc, char ** argv)
{
    seqan3::argument_parser myparser{"build_ibf", argc, argv, false};
    cmd_arguments args{};
    initialize_argument_parser(myparser, args);
    try
    {
         myparser.parse();
    }
    catch (seqan3::argument_parser_error const & ext)
    {
        std::cout << "[Error] " << ext.what() << "\n";
        return -1;
    }

    if (args.gz && args.bz2)
        throw seqan3::argument_parser_error{"Files cannot be both gz and bz2 compressed."};

    args.size.erase(std::remove(args.size.begin(), args.size.end(), ' '), args.size.end());

    if (args.bits != 4096u && !args.size.empty()) // Probably not default. https://github.com/seqan/seqan3/pull/1859
        throw seqan3::argument_parser_error{"Either set --bits or --size."};

    if (!args.size.empty())
    {
        size_t multiplier{};

        switch (std::tolower(args.size.back()))
        {
            case 't':
                multiplier = 8ull * 1024ull * 1024ull * 1024ull * 1024ull;
                break;
            case 'g':
                multiplier = 8ull * 1024ull * 1024ull * 1024ull;
                break;
            case 'm':
                multiplier = 8ull * 1024ull * 1024ull;
                break;
            case 'k':
                multiplier = 8ull * 1024ull;
                break;
            default:
                throw seqan3::argument_parser_error{"Use {k, m, g, t} to pass size. E.g., --size 8g."};
        }

        size_t size{};
        std::from_chars(args.size.data(), args.size.data() + args.size.size() - 1, size);
        size *= multiplier;
        args.bits = size / args.bins;
    }

    run_program(args);
    return 0;
}