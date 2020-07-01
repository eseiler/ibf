#include <mutex>

#include <seqan3/core/debug_stream.hpp>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/range/views/async_input_buffer.hpp>
#include <seqan3/range/views/zip.hpp>
#include <seqan3/search/dream_index/technical_binning_directory.hpp>

struct cmd_arguments
{
    std::filesystem::path query_file{};
    std::filesystem::path ibf_file{};
    std::filesystem::path out_file{"search.out"};
    uint8_t kmer_size{20};
    uint8_t errors{0};
    uint8_t threads{1};
    uint64_t pattern_size{};
    uint8_t parts{1u};
};

struct my_traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

class sync_out
{
public:
    sync_out() = default;
    sync_out(sync_out const &) = default;
    sync_out & operator=(sync_out const &) = default;
    sync_out(sync_out &&) = default;
    sync_out & operator=(sync_out &&) = default;
    ~sync_out() = default;

    sync_out(std::filesystem::path const & path) : file(std::ofstream{path}) {}

    void write(std::string const & data)
    {
        std::lock_guard<std::mutex> lock(write_mutex);
        file << data;
    }

private:
    std::ofstream file;
    std::mutex write_mutex;
};

template <typename t>
void load_ibf(t & tbd, cmd_arguments const & args, size_t const part)
{
    std::filesystem::path ibf_file{args.ibf_file};
    ibf_file += "_" + std::to_string(part);
    std::ifstream is{ibf_file, std::ios::binary};
    cereal::BinaryInputArchive iarchive{is};
    iarchive(tbd);
}

template <typename t>
inline void do_parallel(t && worker, size_t const threads)
{
    std::vector<decltype(std::async(std::launch::async, worker))> tasks;

    for (size_t i = 0; i < threads; ++i)
        tasks.emplace_back(std::async(std::launch::async, worker));

    for (auto && task : tasks)
        task.wait();
}

void run_program_multiple(cmd_arguments const & args)
{
    using std::get;

    seqan3::ibf_config cfg{seqan3::bin_count{64u},
                           seqan3::bin_size{1024u},
                           seqan3::hash_function_count{2u}};

    seqan3::technical_binning_directory tbd{std::vector<std::vector<seqan3::dna4>>{},
                                            seqan3::views::kmer_hash(seqan3::ungapped{args.kmer_size}),
                                            cfg};

    seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::id, seqan3::field::seq>> fin{args.query_file};
    using record_type = typename decltype(fin)::record_type;
    std::vector<record_type> records{};

    for (auto && chunked_records : fin | seqan3::views::chunk(10'000'000))
    {
        records.clear();
        std::ranges::copy(chunked_records, std::cpp20::back_inserter(records));

        load_ibf(tbd, args, 0);

        std::vector<seqan3::counting_vector<uint16_t>> counts(records.size(),
                                                              seqan3::counting_vector<uint16_t>(tbd.bin_count(), 0));

        auto zipped_input = seqan3::views::zip(records, std::views::iota(0u)) | seqan3::views::async_input_buffer(std::max<size_t>(args.threads * 4, 8));

        auto count_task = [&]()
        {
            auto counter = tbd.counting_agent<uint16_t>();

            for (auto && x : zipped_input)
            {
                auto && [id, seq] = get<0>(x);
                size_t i = get<1>(x);
                auto & result = counter.count_query(seq);
                counts[i] += result;
            }
        };

        do_parallel(count_task, args.threads);

        for (size_t const part : std::views::iota(1u, static_cast<unsigned int>(args.parts - 1)))
        {
            load_ibf(tbd, args, part);

            zipped_input = seqan3::views::zip(records, std::views::iota(0u)) | seqan3::views::async_input_buffer(std::max<size_t>(args.threads * 4, 8));
            do_parallel(count_task, args.threads);
        }

        load_ibf(tbd, args, args.parts - 1);
        zipped_input = seqan3::views::zip(records, std::views::iota(0u)) | seqan3::views::async_input_buffer(std::max<size_t>(args.threads * 4, 8));
        sync_out synced_out{args.out_file};
        uint16_t const threshold = args.kmer_size * (args.errors + 1) < args.pattern_size ?
                                       args.pattern_size - (args.kmer_size * (args.errors + 1)) + 1 : 1;

        auto output_task = [&]()
        {
            auto counter = tbd.counting_agent<uint16_t>();
            std::string result_string{};

            for (auto && x : zipped_input)
            {
                auto && [id, seq] = get<0>(x);
                size_t i = get<1>(x);
                auto & result = counter.count_query(seq);
                counts[i] += result;

                result_string.clear();
                result_string += id;
                result_string += '\t';

                size_t current_bin{0};
                for (auto const & count : counts[i])
                {
                    if (count >= threshold)
                    {
                        result_string += std::to_string(current_bin);
                        result_string += ',';
                    }
                    ++current_bin;
                }
                result_string += '\n';
                synced_out.write(result_string);
            }
        };

        do_parallel(output_task, args.threads);
    }
}

void run_program_single(cmd_arguments const & args)
{
    seqan3::ibf_config cfg{seqan3::bin_count{64u},
                           seqan3::bin_size{1024u},
                           seqan3::hash_function_count{2u}};

    seqan3::technical_binning_directory tbd{std::vector<std::vector<seqan3::dna4>>{},
                                            seqan3::views::kmer_hash(seqan3::ungapped{args.kmer_size}),
                                            cfg};

    std::ifstream is{args.ibf_file, std::ios::binary};
    cereal::BinaryInputArchive iarchive{is};
    iarchive(tbd);

    seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::id, seqan3::field::seq>> fin{args.query_file};

    // create the async buffer around the input file
    // spawns a background thread that tries to keep eight records in the buffer
    auto sequence_input_buffer = fin | seqan3::views::async_input_buffer(std::max<size_t>(args.threads * 4, 8));

    sync_out synced_out{args.out_file};
    size_t const threshold = args.kmer_size * (args.errors + 1) < args.pattern_size ? args.pattern_size - (args.kmer_size * (args.errors + 1)) + 1 : 1;

    // create a lambda function that iterates over the async buffer when called
    // (the buffer gets dynamically refilled as soon as possible)
    auto worker = [&] ()
    {
        auto counter = tbd.counting_agent();
        std::string result_string{};

        for (auto && [id, seq] : sequence_input_buffer)
        {
            result_string.clear();
            result_string += id;
            result_string += '\t';

            auto & result = counter.count_query(seq);
            size_t current_bin{0};
            for (auto const & count : result)
            {
                if (count >= threshold)
                {
                    result_string += std::to_string(current_bin);
                    result_string += ',';
                }
                ++current_bin;
            }
            result_string += '\n';
            synced_out.write(result_string);
        }
    };

    do_parallel(worker, args.threads);
}

void initialize_argument_parser(seqan3::argument_parser & parser, cmd_arguments & args)
{
    parser.info.author = "Enrico Seiler";
    parser.info.author = "enrico.seiler@fu-berlin.de";
    parser.info.short_description = "Search reads in a minimizer IBF.";
    parser.info.version = "1.0.0";
    parser.add_positional_option(args.query_file, "Please provide a path the FASTQ file.");
    parser.add_positional_option(args.ibf_file, "Please provide a valid path to an IBF. Parts: Without suffix _0");
    parser.add_option(args.out_file, '\0', "output", "Please provide a valid path to the output.",
                      seqan3::option_spec::DEFAULT);
    parser.add_option(args.kmer_size, '\0', "kmer", "Choose the kmer size.", seqan3::option_spec::DEFAULT,
                      seqan3::arithmetic_range_validator{1, 32});
    parser.add_option(args.threads, '\0', "threads", "Choose the number of threads.", seqan3::option_spec::DEFAULT,
                      seqan3::arithmetic_range_validator{1, 2048});
    parser.add_option(args.parts, '\0', "parts", "Splits the IBF in this many parts. Must be a power of 2.");
    parser.add_option(args.errors, '\0', "error", "Choose the number of errors.", seqan3::option_spec::DEFAULT,
                      seqan3::arithmetic_range_validator{0, 5});
    parser.add_option(args.pattern_size, '\0', "pattern",
                      "Choose the pattern size. Default: Use median of sequence lengths in query file.");
}

int main(int argc, char ** argv)
{
    seqan3::argument_parser myparser{"search_min_ibf", argc, argv, false};
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

    if (!args.pattern_size)
    {
        std::vector<uint64_t> sequence_lengths{};
        seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq>> query_in{args.query_file};
        for (auto & [seq] : query_in | seqan3::views::async_input_buffer(16))
        {
            sequence_lengths.push_back(std::ranges::size(seq));
        }
        std::sort(sequence_lengths.begin(), sequence_lengths.end());
        args.pattern_size = sequence_lengths[sequence_lengths.size()/2];
    }

    if (args.parts == 1)
        run_program_single(args);
    else
        run_program_multiple(args);

    return 0;
}
