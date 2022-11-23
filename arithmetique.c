#include <stdio.h>

long add(unsigned long a[2], unsigned long b[2])
{
    unsigned long c[3];
    int retenue; 
    if ( a[0] + b[0] < a[0])
    // on est en overflow
    {

    }

}

int main(void)
{
    unsigned long a[2], b[2];

    add(a, b);

    return 0;
}