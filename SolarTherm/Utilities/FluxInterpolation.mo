within SolarTherm.Utilities;
function FluxInterpolation
	extends Modelica.Icons.Function;
	input Real nu1;
	input Real nu2;
	input Real nu3;
	input Real nu4;
	input Real nu5;
	input Real ele;
	input Real dni;
	input Real ele_min;
	output Real result;

protected
	Real cons1 = 0.56;
	Real cons2 = 0.87;
	Real cons3 = 1.00;
	Real cons4 = 1.20;
	Real cons5 = 1.39;
	Real x = 1;
	Real dni_clear = 1;
	
algorithm
	if ele > ele_min then
		dni_clear:=1363*0.7^((1./cos(0.5*pi-ele))^0.678);
		x:=dni/dni_clear;
		if x >= cons5 then
			result := nu5;
		elseif x >= cons4 and x<cons5 then
			result := nu4 + (x-cons4)/(cons5-cons4)*(nu5-nu4);
		elseif x >= cons3 and x<cons4 then
			result := nu3 + (x-cons3)/(cons4-cons3)*(nu4-nu3);
		elseif x >= cons2 and x<cons3 then
			result := nu2 + (x-cons2)/(cons3-cons2)*(nu3-nu2);
		elseif x >= cons1 and x<cons2 then
			result := nu1 + (x-cons1)/(cons2-cons1)*(nu2-nu1);
		elseif x < cons1 then
			result := x/cons1*nu1;
		end if;
	else
		dni_clear:=dni;
		x:=0;
		result := 0;
	end if;
end FluxInterpolation;
