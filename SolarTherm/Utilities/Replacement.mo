within SolarTherm.Utilities;
function Replacement
	extends Modelica.Icons.Function;
	input Real r;
	input Integer N;
	input Integer tplant;
	input Real[N] life;
	output Real result;
protected
	Integer res;
	Integer sum;
	Real pvalue;
algorithm
	result:=0;
	sum:=0;
	pvalue:=0;
	for i in 1:N loop
		for j in 1:tplant loop
			if j==1 then
				res:=integer(floor((j+1)/life[i]));
				sum:=0;
				pvalue:=res/(1+r)^j;
			else
				res:=integer(floor((j+1)/life[i])-sum);
				sum:=sum+res;
				pvalue:=pvalue+res/(1+r)^j;
			end if;
		end for;
		result:=result+pvalue;
	end for;
end Replacement;
