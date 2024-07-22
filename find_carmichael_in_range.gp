MILLER_RABIN_K = 40


is_prime_miller_rabin(n, k=MILLER_RABIN_K) = {

	if(n%2 == 0, return(0));
	s = valuation(n-1, 2);
	d = (n-1) / (2^s);
	
	for(i=1, k,
		a = 2  + random(n-4);
		x = Mod((Mod(a, n))^d, n
	);
	
	for(j=1, s,
		y = x^2 % n;
		if(y == 1 && x != 1 && x != n-1, return(0)); 
		x = y
	);

	if(y != 1, return(0)));
	return(1);
}


mark_primes(high) = {

	mprime = List();
	ck = vector(high+1, unused, 1);
	l = floor(sqrt(high));
	
	for(i=2, l,
		if(ck[i]==1,
		  forstep(j=i*i, l+1, i, ck[j]=0;)
	  	   );
		);
		
	for(k=2, l,
		if(ck[k]==1,
			listput(mprime, k);
			);
		);
		
	return(mprime);

}


segmented_sieve(low, high) = {

	primes_output = List();
	if(high <= low, return(primes_output));
	m_primes = mark_primes(high);
	s_primes = vector(high-low+1, unused, 1);


	for(j=1, length(m_primes),
		 i = m_primes[j];
		 lower=floor(low/i);
		 if(lower <= 1,
		 	lower = i+i,
		 	(low % i ) != 0,
		 	lower = (lower*i) + i,
		 	lower=lower*i
		 );
		forstep(k=lower, high, i,
				s_primes[k-low+1]=0);
	);
	
	for(k=low, high,
		if(s_primes[k-low+1]==1,
			listput(primes_output, k)
			);
		);

	return(primes_output);
}


k_th_prime(k) = {

	cache = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29];
	number_of_primes = 10;
	if(k <= number_of_primes, return(cache[k]));

	forstep(i=30, 10^10, 100,
		    sp = Vec(segmented_sieve(i, i+100));
		    number_of_primes += length(sp);
		    if(number_of_primes >= k, break);
	);

	return(sp[k-number_of_primes+length(sp)]);

}


find_carmichaels_in_range(A, B) = {

	(f(m, lam, low, k) = 

		if(k==1,
			low = max(low, floor(A / m) + if(A % m == 0, 0, 1));
			high = min(floor(B / m) + 1, max_p);
			u = (1/m) % lam;
			while(u < low, u +=lam;);
			forstep(p=u, high, lam, if(((((m*p) - 1) % (p-1)  == 0) && (is_prime_miller_rabin(p) == 1 ) ), listput(list, m * p);));
			,
			high = floor(sqrtn(floor(B/m), k) + 1);
			foreach(segmented_sieve(low, high), p, if(gcd(m, p-1)==1, f(m * p, lcm(lam, p-1), p+2, k-1)));
			);
		
	);

	max_p = floor(sqrt(B));
	list = List();
	k=3;
	l = 3 * 5 * 7;
	while(l < B,
		f(1,1,3, k);
		k += 1;
		l *= k_th_prime(k);
		);
	   
	   vecsort(Vec(list));
};



benchmark(A, B, run_n_times) = {

time_list = List();

for(i=1, run_n_times,
	p_0 = getwalltime();
	find_carmichaels_in_range(A, B);
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



