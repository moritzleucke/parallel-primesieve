// compile with:
// g++-7 prime_sieve.cpp -o prime_sieve -std=c++11 -fopenmp
#include <iostream>
#include <cmath>
#include <tuple>
#include <chrono>
#include <omp.h>
using namespace std::chrono;

bool *get_sieve(long long int upper_limit){
    // sieving primes until upper limit, return bool array.
    int sieve_limit = int(sqrt(upper_limit));

    bool *is_prime = new bool[upper_limit];
    // initializing prime array 
    for (long long int i = 0; i < upper_limit; i++) is_prime[i] = true;
    // 1 is not prime
    is_prime[0] = false;

    for (long long int i = 2; i <= sieve_limit; i++) {
        if (is_prime[i-1] == true) {
            for (long long int j = 2*i; j <= upper_limit; j+=i) {
                is_prime[j-1] = false;
            }
        }
    }
    return is_prime;
}

int count_primes(bool *sieve, long long int upper_limit) {
    // count the true values in bool array.
    int n_primes = 0;
    for (long long int i = 0; i < upper_limit; i++) {
        if (sieve[i] == true) {
            n_primes = n_primes + 1;
        }
    }
    return n_primes;
}

std::tuple<int*, int> get_primitive_primes(long long int upper_limit){
    // getting the primitive primes a s array + the number of primitives.
    long long int  primitive_upper_limit = (long long int)sqrt(upper_limit);
    bool *is_prime = get_sieve(primitive_upper_limit);
    int n_primitives = count_primes(is_prime, primitive_upper_limit);
    int *primitive_primes = new int[n_primitives];
    int j = 0;
    for (int i = 0; i < primitive_upper_limit; i++){
        if (is_prime[i] == true) {
            primitive_primes[j] = i+1;
            j++;
        }
    }
    delete[] is_prime;
    return std::make_tuple(primitive_primes, n_primitives);
}

bool *accelerated_sieve(long long int upper_limit, int n_primitives, int *primitive_primes) {
    // accelerated prime sieve until upper_limit using the primitive primes. 
    bool *is_prime = new bool[upper_limit];
    int sieve_limit = int(std::ceil(sqrt(upper_limit)));
    // initialization 
    for (int i = 0; i < sieve_limit; i++) 
        is_prime[i] = false;
    for (int i = 0; i < n_primitives; i++) 
        is_prime[primitive_primes[i]-1] = true;
    #pragma omp parallel for 
    for (long long int i = sieve_limit; i < upper_limit; i++) 
        is_prime[i] = true;
    // sieving
    // TODO: change to balance loops better
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < n_primitives; i++){
        int primitive = primitive_primes[i];
        int multiple_above = int(ceil(sieve_limit/(double)primitive))*primitive;
        for (long long int j = multiple_above; j <= upper_limit; j += primitive){
            is_prime[j-1] = false;
        }
    }
    return is_prime;
}

void print_primes(bool *prime_list, long long int upper_bound){
    for (long long int i = 0; i < upper_bound; i++){
        if (prime_list[i] == true){
            std::cout << i+1 << "\n";
        }
    }
    return;
}

void print_list(int *prime_list, int upper_limit){
    for (int i = 0; i < upper_limit; i++){
        std::cout << prime_list[i] << "\n";
    }
    return;
}

/*
int main() {
    long long int upper_limit;
    std::cout << "enter until which number I should calculate primes: ";
    std::cin >> upper_limit;
    std::cout << "\ncalculating primes until " << upper_limit << " ...\n";

    auto start = high_resolution_clock::now();

    bool *is_prime = get_sieve(upper_limit);

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    double du = double(duration.count())/1000000.0;
    std::cout << "execution time: " << du << "\n";

    int n_primes = count_primes(is_prime, upper_limit);
    std::cout << "number of primes until " << upper_limit << " is " << n_primes << ".\n";
    std::cout << sizeof(bool) << ", " << sizeof(char);

    delete[] is_prime;
}
*/

int main(){
    int *primitive_primes;
    bool *is_prime;
    int n_primitives;
    long long int upper_limit;
    std::cout << "enter until which number I should calculate primes: ";
    std::cin >> upper_limit;
    std::cout << "\ncalculating primes until " << upper_limit << " ...\n";


    auto start = high_resolution_clock::now();
    std::tie(primitive_primes, n_primitives) = get_primitive_primes(upper_limit);
    is_prime = accelerated_sieve(upper_limit, n_primitives, primitive_primes);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    double du = double(duration.count())/1000000.0;
    std::cout << "execution time: " << du << "\n";
    #pragma omp parallel
    {
        #pragma omp single
        std::cout << "executed on " << omp_get_num_threads() << " threads\n";
    }
    
    
    int n_primes = count_primes(is_prime, upper_limit);
    std::cout << "number of primes until " << upper_limit << " is " << n_primes << ".\n";
    std::cout << "sieve limit: " << int(sqrt(upper_limit)) << "\n";
    std::cout << "number of scanned primes: " << n_primitives << "\n";

    
    // print_list(primitive_primes, n_primitives);
    // std::cout << "\n\n";
    // print_primes(is_prime, upper_limit);

    delete[] primitive_primes;
    delete[] is_prime;
}
