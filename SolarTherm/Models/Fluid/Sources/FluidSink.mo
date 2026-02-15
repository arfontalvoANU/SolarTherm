within SolarTherm.Models.Fluid.Sources;
model FluidSink "Infinite fluid sink"
    replaceable package Medium = Modelica.Media.Interfaces.PartialMedium;

    Modelica.Fluid.Interfaces.FluidPort_a port_a(
        redeclare package Medium=Medium,
        m_flow(min=0))
        annotation (
            Placement(
                transformation(
                    extent={{-110,-12},{-90,8}},
                    rotation=0),
                iconTransformation(
                    extent={{-110,-12},{-90,8}})));
equation
    port_a.h_outflow = 0; // shouldn't flow backwards anyway
    port_a.p=101325;
    port_a.Xi_outflow=inStream(port_a.Xi_outflow);

annotation(
    Icon(graphics={
        Ellipse(
            extent={{-100,100},{100,-100}},
            fillPattern=FillPattern.Sphere,
            fillColor={0,127,255}),
        Text(
            extent={{-150,110},{150,150}},
            textString="%name",
            lineColor={0,0,255})
        }),
    Documentation(info =
        "<html>
        <ul>
        <li> <i>Dec 2020</i> by A. Shirazi:<br> Resealed first version. </li>
        <li> <i>Feb 2026</i> by A. Fontalvo:<br> Icon improvement. </li>
        </ul>
        </html>"));end FluidSink;
