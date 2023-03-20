#include <iostream>
#include <iomanip>
#include <chrono>
#include <cstring>
#include <cstdio>
#include <unistd.h>
#include <sys/ioctl.h>

namespace bcolors {
    const std::string BLACK = "\033[0;30m";
    const std::string RED = "\033[0;31m";
    const std::string GREEN = "\033[0;32m";
    const std::string BROWN = "\033[0;33m";
    const std::string BLUE = "\033[0;34m";
    const std::string PURPLE = "\033[0;35m";
    const std::string CYAN = "\033[0;36m";
    const std::string LIGHT_GRAY = "\033[0;37m";
    const std::string DARK_GRAY = "\033[1;30m";
    const std::string LIGHT_RED = "\033[1;31m";
    const std::string LIGHT_GREEN = "\033[1;32m";
    const std::string YELLOW = "\033[1;33m";
    const std::string LIGHT_BLUE = "\033[1;34m";
    const std::string LIGHT_PURPLE = "\033[1;35m";
    const std::string LIGHT_CYAN = "\033[1;36m";
    const std::string LIGHT_WHITE = "\033[1;37m";
    const std::string BOLD = "\033[1m";
    const std::string FAINT = "\033[2m";
    const std::string ITALIC = "\033[3m";
    const std::string UNDERLINE = "\033[4m";
    const std::string BLINK = "\033[5m";
    const std::string NEGATIVE = "\033[7m";
    const std::string CROSSED = "\033[9m";
    const std::string ENDC = "\033[0m";
}

class ScreenOutput {
    private:
        int n_columns;
        std::chrono::time_point<std::chrono::high_resolution_clock> T;

    public:
        ScreenOutput() {
            struct winsize size;
            ioctl(STDOUT_FILENO, TIOCGWINSZ, &size);
            n_columns = std::min(size.ws_col, static_cast<short unsigned int>(70));
            T = std::chrono::high_resolution_clock::now();
        }

        void print_header() {
            std::string line(n_columns, '=');
            std::string title = "Gamma-SMC";
            std::cout << bcolors::BLUE << bcolors::BOLD << line << std::endl;
            std::cout << std::setw((n_columns+title.length())/2) << title << std::endl;
            std::cout << line << bcolors::ENDC << std::endl;
        }

        void print_title(std::string text) {
            T = std::chrono::high_resolution_clock::now();
            std::string line(n_columns, '-');
            std::cout << bcolors::BLUE << line << std::endl;
            std::cout << " " << text << std::endl;
            std::cout << line << bcolors::ENDC << std::endl;
        }

        void print_subtitle(std::string text) {
            T = std::chrono::high_resolution_clock::now();
            std::cout << bcolors::BLUE << "--- " << text << bcolors::ENDC << std::endl;
        }

        void print_item(std::string text) {
            std::cout << "  * " << text << std::endl;
        }

        void print_done(double secs = -1) {
            if (secs < 0) {
                auto now = std::chrono::high_resolution_clock::now();
                secs = std::chrono::duration_cast<std::chrono::duration<double>>(now - T).count();
            }
            std::cout << std::string(n_columns - 17 - std::to_string(secs).length(), ' ') << bcolors::GREEN << "..done (" << std::fixed << std::setprecision(2) << secs << " secs)." << bcolors::ENDC  << std::endl;
        }
};