R_truncation_rule_domaining_single=function(sd_vector)
{
  if(sd_vector[1]<0){domain=1}
  else if (sd_vector[2]<0){domain=2}
  else if (sd_vector[3]<0){domain=3}else{domain=4}
  return(domain)
}

R_truncation_rule_domaining=function(sd_matrix)
{
  return(apply(sd_matrix,1,R_truncation_rule_domaining_single))
}



