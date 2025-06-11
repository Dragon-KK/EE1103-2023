#include <stdio.h>
/**
 * Calculates x^y vro
 * NOTE: y must be positive boss
 * 
 * @param x Base
 * @param y Exponent
*/
long long int powerBro(long long int x, unsigned long long int y){
    if (y == 0) return 1;
    else if (y == 1) return x;

    if (y%2 == 0){
        long long int tmp = powerBro(x, y/2);
        return tmp*tmp;
    }
    else{
        long long int tmp = powerBro(x, y/2);

        return tmp*tmp*x;
    }
}


int main(){
    printf("%Ld", powerBro(5, 5));
}