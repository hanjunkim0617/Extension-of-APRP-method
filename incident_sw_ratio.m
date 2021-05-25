function QS=incident_sw_ratio(clt,a_clr,a_oc,mu_clr,mu_cld,ga_clr,ga_cld,flag_model,ar,rr)

mu_oc=-(1-mu_clr).*(1-mu_cld)+1;
ga_oc=-(1-ga_clr).*(1-ga_cld)+1;

if flag_model == 1
	QS_clr = (1-mu_clr).*(1-ga_clr) ./ (1-(1-ar.*mu_clr).*a_clr.*rr.*ga_clr);
	QS_oc = (1-mu_oc).*(1-ga_oc) ./ (1-(1-ar.*mu_oc).*a_oc.*rr.*ga_oc);
elseif flag_model == 2
	QS_clr = (1-mu_clr-ga_clr) ./ (1-a_clr.*rr.*ga_clr);
	QS_oc = (1-mu_oc-ga_oc) ./ (1-a_oc.*rr.*ga_oc);
end

QS=(1-clt).*QS_clr+clt.*QS_oc;

end

