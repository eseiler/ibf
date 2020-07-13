#include <chrono>
#include <mutex>

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
    bool write_time{false};
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
void load_ibf(t & tbd, cmd_arguments const & args, size_t const part, double & ibf_io_time)
{
    std::filesystem::path ibf_file{args.ibf_file};
    ibf_file += "_" + std::to_string(part);
    std::ifstream is{ibf_file, std::ios::binary};
    cereal::BinaryInputArchive iarchive{is};
    auto start = std::chrono::high_resolution_clock::now();
    iarchive(tbd);
    auto end = std::chrono::high_resolution_clock::now();
    ibf_io_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
}

template <typename t>
inline void do_parallel(t && worker, size_t const num_records, size_t const threads, double & compute_time)
{
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<decltype(std::async(std::launch::async, worker, size_t{}, size_t{}))> tasks;
    size_t const records_per_thread = num_records / threads;

    for (size_t i = 0; i < threads; ++i)
    {
        size_t const start = records_per_thread * i;
        size_t const end = i == (threads-1) ? num_records: records_per_thread * (i+1);
        tasks.emplace_back(std::async(std::launch::async, worker, start, end));
    }

    for (auto && task : tasks)
        task.wait();

    auto end = std::chrono::high_resolution_clock::now();
    compute_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
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

    double ibf_io_time{0.0};
    double reads_io_time{0.0};
    double compute_time{0.0};

    auto cereal_worker = [&] ()
    {
        load_ibf(tbd, args, 0, ibf_io_time);
    };

    for (auto && chunked_records : fin | seqan3::views::chunk((1ULL<<20)*10))
    {
        auto cereal_handle = std::async(std::launch::async, cereal_worker);

        records.clear();
        auto start = std::chrono::high_resolution_clock::now();
        std::ranges::move(chunked_records, std::cpp20::back_inserter(records));
        auto end = std::chrono::high_resolution_clock::now();
        reads_io_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();

        std::vector<seqan3::counting_vector<uint16_t>> counts(records.size(),
                                                              seqan3::counting_vector<uint16_t>(tbd.bin_count(), 0));

        auto count_task = [&](size_t const start, size_t const end)
        {
            auto counter = tbd.counting_agent<uint16_t>();
            size_t counter_id = start;

            for (auto && [id, seq] : records | seqan3::views::slice(start, end))
            {
                auto & result = counter.count_query(seq);
                counts[counter_id++] += result;
            }
        };

        cereal_handle.wait();
        do_parallel(count_task, records.size(), args.threads, compute_time);

        for (size_t const part : std::views::iota(1u, static_cast<unsigned int>(args.parts - 1)))
        {
            load_ibf(tbd, args, part, ibf_io_time);
            do_parallel(count_task, records.size(), args.threads, compute_time);
        }

        load_ibf(tbd, args, args.parts - 1, ibf_io_time);
        sync_out synced_out{args.out_file};
        uint16_t const threshold = args.kmer_size * (args.errors + 1) < args.pattern_size ?
                                       args.pattern_size - (args.kmer_size * (args.errors + 1)) + 1 : 1;

        auto output_task = [&](size_t const start, size_t const end)
        {
            auto counter = tbd.counting_agent<uint16_t>();
            size_t counter_id = start;
            std::string result_string{};

            for (auto && [id, seq] : records | seqan3::views::slice(start, end))
            {
                auto & result = counter.count_query(seq);
                counts[counter_id] += result;

                result_string.clear();
                result_string += id;
                result_string += '\t';

                size_t current_bin{0};
                for (auto && count : counts[counter_id++])
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

        do_parallel(output_task, records.size(), args.threads, compute_time);
    }

    if (args.write_time)
    {
        std::filesystem::path file_path{args.out_file};
        file_path += ".time";
        std::ofstream file_handle{file_path};
        file_handle << "IBF I/O\tReads I/O\tCompute\n";
        file_handle << std::fixed
                    << std::setprecision(2)
                    << ibf_io_time << '\t'
                    << reads_io_time << '\t'
                    << compute_time;
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

    double ibf_io_time{0.0};
    double reads_io_time{0.0};
    double compute_time{0.0};

    auto cereal_worker = [&] ()
    {
        auto start = std::chrono::high_resolution_clock::now();
        iarchive(tbd);
        auto end = std::chrono::high_resolution_clock::now();
        ibf_io_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
    };
    auto cereal_handle = std::async(std::launch::async, cereal_worker);

    seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::id, seqan3::field::seq>> fin{args.query_file};
    using record_type = typename decltype(fin)::record_type;
    std::vector<record_type> records{};

    sync_out synced_out{args.out_file};
    size_t const threshold = args.kmer_size * (args.errors + 1) < args.pattern_size ?
                                args.pattern_size - (args.kmer_size * (args.errors + 1)) + 1 : 1;

    auto worker = [&] (size_t const start, size_t const end)
    {
        auto counter = tbd.counting_agent<uint16_t>();
        std::string result_string{};

        for (auto && [id, seq] : records | seqan3::views::slice(start, end))
        {
            result_string.clear();
            result_string += id;
            result_string += '\t';

            auto & result = counter.count_query(seq);
            size_t current_bin{0};
            for (auto && count : result)
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

    for (auto && chunked_records : fin | seqan3::views::chunk((1ULL<<20)*10))
    {
        records.clear();
        auto start = std::chrono::high_resolution_clock::now();
        std::ranges::move(chunked_records, std::cpp20::back_inserter(records));
        auto end = std::chrono::high_resolution_clock::now();
        reads_io_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();

        cereal_handle.wait();

        do_parallel(worker, records.size(), args.threads, compute_time);
    }

    if (args.write_time)
    {
        std::filesystem::path file_path{args.out_file};
        file_path += ".time";
        std::ofstream file_handle{file_path};
        file_handle << "IBF I/O\tReads I/O\tCompute\n";
        file_handle << std::fixed
                    << std::setprecision(2)
                    << ibf_io_time << '\t'
                    << reads_io_time << '\t'
                    << compute_time;
    }
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
    parser.add_flag(args.write_time, '\0', "time", "Write timing file.", seqan3::option_spec::ADVANCED);
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
