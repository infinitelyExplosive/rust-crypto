
from random import randint
from math import gcd


def find_func(a):
    b = randint(2,60)
    for c in range(30,50):
        for d in range(150,250):
            for x in range(1,40):
                for y in range(1,60):
                    if a + x * b + y * c + x * y * d == 0:
                        if gcd(a, gcd(b, gcd(c, d))) == 1:
                            cap_w = max([a, b*(x+1), c*(y+1), d*(x+1)*(y+1)])
                            if (x+1)*(y+1) < cap_w**(2/3 - 1/6):
                                msg = f'{a},{b},{c},{d},{x},{y}'
                                print(msg)
                                return msg
    return ''

with open('polys.txt', 'w') as f:
    # f.write('a,b,c,d,x,y\n')
    a = -50000
    while a < -20000:
        result = find_func(a)
        if len(result) > 0:
            f.write(result + '\n')
        a += randint(1,100)