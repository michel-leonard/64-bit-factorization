This  integer factorization algorithm is released into the **public domain**, allowing for unrestricted use, modification, and distribution.

# 64-bit Integer Factorization Algorithm

Your C implementation of a 64-bit integer factorization algorithm utilizes Pollard's Rho method for factoring large numbers.

### Key Features

- **C99 Compliant**: The code is written in pure C99 with no external dependencies, ensuring compatibility and ease of integration into other projects.
- **Pollard's Rho**: For larger numbers, the algorithm leverages Pollard's Rho method to efficiently find factors.
- **No Dependencies**: The implementation is completely standalone and does not require any external libraries or includes, making it lightweight and easy to use.

### Getting Started

To use the factorization algorithm, simply include the provided source code in your project and call the main factorization function with the desired 64-bit integer. Compile using `gcc -Wall -Werror -Wextra -pedantic -O3 -std=c99 factor.c -o fac.exe`.
