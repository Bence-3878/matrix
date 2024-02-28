#include <iostream>
#include "Vektor.h"
#include "Matrix.h"

int main() {
    /* 3 * x - 4 * y + 0 * z = 6
     * 1 * x + 1 * y + 4 * z = 5
     * 0 * x - 1 * y + 5 * z = -4
     */
    std::cout << "hány egyenlet van?" << std::endl;
    int n,m;
    std::cin >> n;
    std::cout << "hány ismeretlen van?" << std::endl;
    std::cin >> m;
    Matrix<double> egyenlet(n,m);
    char* ism = new char[m];
    for (int i = 0; i < n; ++i) {
        std::cin >> egyenlet[]
    }
    return 0;
}
