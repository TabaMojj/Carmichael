find_carmichaels_in_range_builtin_functions(A, B) = {

max_p = floor(sqrt(B));
list = List();
(f(m, lam, low, k) = 

if(k==1,
	low = max(low, floor(A / m) + if(A % m == 0, 0, 1));
	high = min(floor(B / m) + 1, max_p);
	u = (1/m) % lam;
	while(u < low, u +=lam;);
	forstep(p=u, high, lam,
	        if(((((m*p) - 1) % (p-1)  == 0) && (isprime(p)) ), listput(list, m * p);));,
	high = sqrtn(floor(B/m), k) + 1;
	forprime(p=low, high, if(gcd(m, p-1)==1, f(m * p, lcm(lam, p-1), p+2, k-1)));
	);
   );
   
   k=3;
   l = 3 * 5 * 7;
   while(l < B,
    f(1,1,3, k);
     k += 1;
     l *= prime(k);
     
      );
   
   vecsort(Vec(list));
   };
 
 
benchmark(A, B, run_n_times) = {

time_list = List();

for(i=1, run_n_times,
	p_0 = getwalltime();
	find_carmichaels_in_range_builtin_functions(A, B);
	p_1=getwalltime();
	listput(time_list, (p_1-p_0)*0.001);
	);
	
time_vec = Vec(time_list);

sd=0;
time_mean = vecsum(time_vec)/run_n_times;
foreach(time_vec, t,  sd +=  ((t-time_mean)^2/(run_n_times-1));  );
sd = sqrt(sd);
printf("Mean: %f", time_mean );
print("\n");
printf("Standard Deviation: %f", sd );
}
