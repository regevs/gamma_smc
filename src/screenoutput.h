#include <iostream>
#include <iomanip>
#include <chrono>
#include <cstring>
#include <cstdio>
#include <unistd.h>
#include <sys/ioctl.h>

namespace bcolors {
    const string BLACK = "\033[0;30m";
    const string RED = "\033[0;31m";
    const string GREEN = "\033[0;32m";
    const string BROWN = "\033[0;33m";
    const string BLUE = "\033[0;34m";
    const string PURPLE = "\033[0;35m";
    const string CYAN = "\033[0;36m";
    const string LIGHT_GRAY = "\033[0;37m";
    const string DARK_GRAY = "\033[1;30m";
    const string LIGHT_RED = "\033[1;31m";
    const string LIGHT_GREEN = "\033[1;32m";
    const string YELLOW = "\033[1;33m";
    const string LIGHT_BLUE = "\033[1;34m";
    const string LIGHT_PURPLE = "\033[1;35m";
    const string LIGHT_CYAN = "\033[1;36m";
    const string LIGHT_WHITE = "\033[1;37m";
    const string BOLD = "\033[1m";
    const string FAINT = "\033[2m";
    const string ITALIC = "\033[3m";
    const string UNDERLINE = "\033[4m";
    const string BLINK = "\033[5m";
    const string NEGATIVE = "\033[7m";
    const string CROSSED = "\033[9m";
    const string ENDC = "\033[0m";
}

class ScreenOutput {
    private:
        int n_columns;
        std::chrono::time_point<std::chrono::high_resolution_clock> T;

    public:
        ScreenOutput() {
            n_columns = 70; 
            T = std::chrono::high_resolution_clock::now();
        }

        void print_header() {
            string line(n_columns, '=');
            string title = "Gamma-SMC";
            cout << endl;
            cout << bcolors::BLUE << bcolors::BOLD << title << bcolors::ENDC << endl;
            cout << endl;
        }

        void print_title(string text) {
            T = std::chrono::high_resolution_clock::now();
            string line(n_columns, '-');
            cout << bcolors::BLUE << line << endl;
            cout << " " << text << endl;
            cout << line << bcolors::ENDC << endl;
        }

        void print_subtitle(string text) {
            T = std::chrono::high_resolution_clock::now();
            cout << bcolors::BLUE << "--- " << text << bcolors::ENDC << endl;
        }

        void print_item(string text) {
            cout << "  * " << text << endl;
        }

        void print_done(double secs = -1) {
            if (secs < 0) {
                auto now = std::chrono::high_resolution_clock::now();
                secs = std::chrono::duration_cast<std::chrono::duration<double>>(now - T).count();
            }
            cout << string(n_columns - 17 - std::to_string(secs).length(), ' ') << bcolors::GREEN << "..done (" << std::fixed << std::setprecision(2) << secs << " secs)." << bcolors::ENDC  << endl;
        }
};