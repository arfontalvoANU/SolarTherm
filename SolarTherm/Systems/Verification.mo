within SolarTherm.Systems;
package Verification
extends Modelica.Icons.Package;
  record CaseData
    extends Modelica.Icons.Record;
    import SI = Modelica.SIunits;
    parameter SI.Temperature T_max;
    parameter SI.Temperature T_min;
    parameter SI.Length H_tank;
    parameter SI.Diameter D_tank;
    parameter SI.CoefficientOfHeatTransfer U_loss_tank;
    parameter SI.MassFlowRate m_flow;
    parameter SI.Length d_p;
  end CaseData;

  record ScenarioBank
    extends Modelica.Icons.Record;
    parameter SI.MassFlowRate m_flow_1 = 0.225 * 0.25*Modelica.Constants.pi*0.148^2;
    parameter CaseData case1(T_max=823, T_min=293, H_tank=1.2, D_tank=0.148, U_loss_tank=0.678, m_flow=m_flow_1, d_p=0.02);
    parameter CaseData case2(T_max=1173, T_min=613, H_tank=13.32, D_tank=2.66, U_loss_tank=0, m_flow=13.7, d_p=0.02);
  end ScenarioBank;
end Verification;