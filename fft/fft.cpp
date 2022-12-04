#include <cmath>
#include <vector>
#include <complex>
#include <iostream>
#include <matplot/matplot.h>

typedef std::vector<std::complex<double> > vec_complex;
typedef std::vector<double> vec_real;

std::complex<double> omega(int n) {
    return exp(std::complex<double>(0, -2 * M_PI / n));
}

// vec_real dft(vec_real input) {
//     vec_real result = vec_real();
// }

/**
 * @brief Fast Fourier Transform
 * 
 * @param input Input signal (complex)
 * @return vec_complex 
 */
vec_complex fft(vec_complex input) {
    size_t n = input.size();
    if (n == 1) {
        return input;
    }

    vec_complex p_e = vec_complex((n+1)/2, 0);
    vec_complex p_o = vec_complex(n/2, 0);

    for (int i = 0, k = 0; i < n-1; i += 2, k++) {
        p_e[k] = input[i];
        p_o[k] = input[i+1];
    }

    vec_complex y_e = fft(p_e);
    vec_complex y_o = fft(p_o);

    vec_complex result = vec_complex(n, 0);
    for (int j = 0; j < n/2; j++) {
        result[j]     = y_e[j] + pow(omega(n), j)*y_o[j];
        result[j+n/2] = y_e[j] - pow(omega(n), j)*y_o[j];
    }

    return result;
}

vec_complex to_complex(vec_real in) {
    vec_complex out;
    for (int i = 0; i < in.size(); i++) {
        out.push_back(std::complex<double>(in[i], 0));
    }
    return out;
}

vec_real to_real(vec_complex in) {
    vec_real r = vec_real(in.size(), 0);
    for (int i = 0; i < in.size(); i++) {
        r[i] = in[i].real();
    }
    return r;
}

vec_real fft_real(vec_real input) {
    vec_complex c = fft(to_complex(input));
    vec_real r = vec_real(input.size(), 0);
    for (int i = 0; i < input.size(); i++) {
        r[i] = c[i].real();
    }
    return r;
}

int main() {
    std::cout << "FFT" << std::endl;

    // test 1
    vec_real data {1,2,3,4,4,3,2,1};
    vec_complex res1 = fft(to_complex(data));
    std::cout << "-- result -- " << std::endl;
    for (int i = 0; i < res1.size(); i++) {
        std::cout << res1[i] << std::endl;
    }

    // test 2
    vec_real data2;
    for (int i = 0; i < 1000; i++) {
        data2.push_back(20*sin(i/1000));
    }
    matplot::plot(data2, fft_real(data2), "x");
    matplot::show();
}
