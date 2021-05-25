function A=albedo(clt,a_clr,a_oc,mu_clr,mu_cld,ga_clr,ga_cld,flag_model,ar,rr)

mu_oc=-(1-mu_clr).*(1-mu_cld)+1;
ga_oc=-(1-ga_clr).*(1-ga_cld)+1;


if flag_model == 1
	A_clr=(1-mu_clr).*ga_clr+(1-mu_clr).*(1-ga_clr).*(1-ar.*mu_clr).*(1-rr.*ga_clr).*a_clr./(1-(1-ar.*mu_clr).*a_clr.*rr.*ga_clr);
	A_oc=(1-mu_oc).*ga_oc+(1-mu_oc).*(1-ga_oc).*(1-ar.*mu_oc).*(1-rr.*ga_oc).*a_oc./(1-(1-ar.*mu_oc).*a_oc.*rr.*ga_oc);
elseif flag_model == 2
	A_clr=ga_clr+(1-mu_clr-ga_clr).*(1-ar.*mu_clr-rr.*ga_clr).*a_clr./(1-a_clr.*rr.*ga_clr);
	A_oc=ga_oc+(1-mu_oc-ga_oc).*(1-ar.*mu_oc-rr.*ga_oc).*a_oc./(1-a_oc.*rr.*ga_oc);
end

A=(1-clt).*A_clr+clt.*A_oc;

end

