within SolarTherm.Models.Fluid.Sources;

model FluidSink2 "Fluid sink but does not fix any pressure or mass-fraction value."
    replaceable package Medium = Modelica.Media.Interfaces.PartialMedium;
    Modelica.Fluid.Interfaces.FluidPort_a port_a(
        redeclare package Medium=Medium,
        m_flow(min=0))
        annotation (Placement(
            visible = true,
            transformation(
                extent={{-110,-12},{-90,8}},
                rotation=0),
            iconTransformation(
                origin = {-100, 8.88178e-16},
                extent = {{-6, -6}, {6, 6}},
                rotation = 0)));
           
equation
    port_a.h_outflow = inStream(port_a.h_outflow); // shouldn't flow backwards anyway

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
        <li> <i>Dec 2020</i> by Z. Kee:<br> Resealed first version. </li>
        <li> <i>Feb 2026</i> by A. Fontalvo:<br> Icon improvement. </li>
        </ul>
        </html>"));
end FluidSink2;