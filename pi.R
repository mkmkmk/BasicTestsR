siz = 10000

min_diff = 1e6
min_a = 0
min_b = 0

a=355
b=113

for (a in 1:siz)
{
    for (b in 1:siz)
    {
        if (abs(a - 3 * b) > b) 
            next
        diff = abs(a / b - pi)
        if (diff<min_diff)
        {
            print(data.frame(a, b, diff))
            min_a = a
            min_b = b
            min_diff =  diff
        }
    }
}


print(data.frame(min_a, min_b, min_diff))

pi
22/7
355/113
3 + 16/113
16/113
16 / (1 + 7 * 16)
3 + 1 / (7 + 1/16)

1 + 15/7
3 + 1/7

355/113 * (1 - 3e-4/3533) - pi
355/113 - pi
sprintf("%.20g" , 355/113)
sprintf("%.20g" , 355/113 * (1 - 3e-4/3533))
sprintf("%.20g" , 312689/99532 )
sprintf("%.60g" , pi)
sprintf("%.60g" , 245850922 / 78256779)

104348/33215 - pi
312689/99532 - pi
833719/265381 - pi
1146408/364913 - pi

245850922 / 78256779 - pi
2549491779/811528438 - pi


