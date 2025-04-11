gsindex<-function(n)	{
	if (sum(n > 2) < 2)
		return(NA)
	length(n) / mean((n - 1) / n)
}
