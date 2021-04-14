#  Fermat's Little Theorem

((1:30)^1) %% 1
((1:30)^2) %% 2
((1:30)^3) %% 3
((1:30)^5) %% 5
((1:30)^7) %% 7
((1:30)^11) %% 11

p = 7
mx = 100
all ((1:mx)^p %% p == 1:mx %% p)

p = 11
mx = 20
all ((1:mx)^p %% p == 1:mx %% p)


mx = 20
for(p in 1:25)
{
    # (ograniczenie wg możliwej dokładności reprezentacji float)
    mx = min(1000, trunc(1e15^(1/p)))
    ferm = all((1:mx)^p %% p == 1:mx %% p) 
    cat(sprintf("p = %2g, F = %5s  (mx = %g) \n", p, ferm, mx))
}

