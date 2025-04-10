function l= lambda_c(lambda_s, Q0, q,  c)
l = (lambda_s+Q0/q)*(c/(c+1)) -Q0/q;
if l<0
    l=0;
end