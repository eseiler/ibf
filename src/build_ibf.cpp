#include <seqan3/std/charconv>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/range/views/async_input_buffer.hpp>
#include <seqan3/range/views/get.hpp>
#include <seqan3/range/views/persist.hpp>
#include <seqan3/range/views/move.hpp>
#include <seqan3/search/dream_index/technical_binning_directory.hpp>

struct my_traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

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
    bool compressed{false};
    bool gz{false};
    bool bz2{false};
    uint8_t parts{1u};
};

template <bool compressed>
struct ibf_builder
{
    cmd_arguments const * arguments;

    template <typename view_t = std::ranges::empty_view<int>>
    auto construct(view_t && restrict_view = std::ranges::empty_view<int>())
    {
        std::string extension{".fasta"};
        if (arguments->gz)
            extension += ".gz";
        if (arguments->bz2)
            extension += ".bz2";

        using sequence_file_t = seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq>>;

        std::vector<sequence_file_t> technical_bins;
        technical_bins.reserve(arguments->bins);

        for (size_t i = 0; i < arguments->bins; ++i)
            technical_bins.emplace_back(arguments->bin_path / ("bin_" + std::to_string(i) + extension));

        seqan3::ibf_config cfg{seqan3::bin_count{arguments->bins},
                               seqan3::bin_size{arguments->bits / arguments->parts},
                               seqan3::hash_function_count{arguments->hash},
                               arguments->threads};

        if constexpr (std::same_as<view_t, std::ranges::empty_view<int>>)
        {
            return seqan3::technical_binning_directory{std::move(technical_bins),
                                                       seqan3::views::kmer_hash(seqan3::ungapped{arguments->k}),
                                                       cfg};
        }
        else
        {
            return seqan3::technical_binning_directory{std::move(technical_bins),
                                                       seqan3::views::kmer_hash(seqan3::ungapped{arguments->k}) |
                                                           restrict_view,
                                                       cfg};
        }
    }

    template <typename view_t = std::ranges::empty_view<int>>
        requires !compressed
    auto ibf(view_t && restrict_view = std::ranges::empty_view<int>())
    {
        return construct(std::move(restrict_view));
    }

    template <typename view_t = std::ranges::empty_view<int>>
        requires compressed
    auto ibf(view_t && restrict_view = std::ranges::empty_view<int>())
    {
        auto tmp = construct(std::move(restrict_view));

        return seqan3::technical_binning_directory<seqan3::data_layout::compressed,
                                                   typename decltype(tmp)::hash_adaptor_t,
                                                   seqan3::dna4>{std::move(tmp)};
    }
};

template <bool compressed>
void run_program(cmd_arguments & args)
{
    ibf_builder<compressed> generator{&args};

    if (args.parts == 1u)
    {
        auto tbd = generator.ibf();
        std::ofstream os{args.out_path, std::ios::binary};
        cereal::BinaryOutputArchive oarchive{os};
        oarchive(tbd);
    }
    else
    {
        std::vector<std::vector<size_t>> association(args.parts);
        size_t next_power_of_four{4u};

        if (args.parts == 4u) // one-to-one
        {
            for (size_t i : std::views::iota(0u, args.parts))
                association[i] = std::vector<size_t>{i};
        }
        else if (args.parts == 2u) // More than 1 prefix per part
        {
            association[0] = std::vector<size_t>{0, 1};
            association[1] = std::vector<size_t>{2, 3};
        }
        else // Multiple prefixes per part
        {
            // How long must the suffix be such that 4^suffix_length >= args.parts
            size_t suffix_length{0};
            for (; 0b100 << (2 * suffix_length) < args.parts; ++suffix_length) {}
            next_power_of_four = 0b100 << (2 * suffix_length);

            size_t const prefixes_per_part = next_power_of_four / args.parts;

            for (size_t i : std::views::iota(0u, next_power_of_four))
                association[i/prefixes_per_part].push_back(i);
        }

        for (size_t part : std::views::iota(0u, args.parts))
        {
            size_t const mask{next_power_of_four - 1};
            auto filter_view = std::views::filter([&] (auto && hash)
                { return std::ranges::find(association[part], hash & mask) != association[part].end(); });

            auto tbd = generator.ibf(filter_view);
            std::filesystem::path out_path{args.out_path};
            out_path += "_" + std::to_string(part);
            std::ofstream os{out_path, std::ios::binary};
            cereal::BinaryOutputArchive oarchive{os};
            oarchive(tbd);
        }
    }
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
    parser.add_option(args.parts, '\0', "parts", "Splits the IBF in this many parts. Must be a power of 2.");
    parser.add_option(args.hash, '\0', "hash", "Choose the number of hashes.", seqan3::option_spec::DEFAULT,
                      seqan3::arithmetic_range_validator{1, 4});
    parser.add_flag(args.compressed, '\0', "compressed", "Build a compressed IBF.");
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
        args.bits = size / (((args.bins + 63) >> 6) << 6);
    }

    if (args.compressed)
        run_program<true>(args);
    else
        run_program<false>(args);
    return 0;
}
