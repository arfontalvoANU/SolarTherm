block TestSTTable
	import SolarTherm.Utilities.Tables.STTable;
	import SolarTherm.Utilities.Tables.tableGrid;
	import SolarTherm.Utilities.Tables.tableCatromInterp;
	STTable table = STTable("resources/tests/sym_test_file1.csv", delim=",");
	parameter Integer res(fixed=false);
	parameter Real val(fixed=false);
initial algorithm
	res := tableGrid(table, 10, 24, "E");
	val := tableCatromInterp(table, 0., 0., 360./10, 360./24);
end TestSTTable;