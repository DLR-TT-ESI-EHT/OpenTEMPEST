within ;
package OpenTEMPEST "Transient Electrochemical Model for Process and Energy SysTems"

  package UsersGuide "User's Guide"
    extends Modelica.Icons.Information;

    class Overview "Overview of the modelling levels"
      extends Modelica.Icons.Information;
      annotation (Documentation(info="<html>
<p>OpenTEMPEST builds large models from small ones. Every level is a model of its own that you can simulate,
replace or reuse in a level above. Depending on the question, each level can be modelled in detail or with a
proven simplification.</p>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"4\">
<tr><th>Level</th><th>What it describes</th><th>Main models</th></tr>
<tr><td>1. Electrochemistry</td>
    <td>Voltage, reaction rates and reaction heat as a function of temperature, gas composition, pressure and
    current density. Lumped (ASR) or detailed (Butler-Volmer, ohmic and diffusion losses).</td>
    <td><a href=\"modelica://OpenTEMPEST.SOC.Electrochem\">SOC.Electrochem</a></td></tr>
<tr><td>2. Cell (repeating unit)</td>
    <td>PEN (positive electrode, electrolyte, negative electrode), fuel and air channels and interconnects with
    mass, momentum and energy balances, discretised with finite volumes in 1D (co-flow) or 2D (cross-flow).</td>
    <td><a href=\"modelica://OpenTEMPEST.SOC.Cell\">SOC.Cell</a></td></tr>
<tr><td>3. Stack</td>
    <td>Cells stacked between end plates with gas distribution and heat transfer to the surroundings, from a
    lumped 0D stack to stacks of 2D cells, with simplified stacks that resolve only a few cells in detail.</td>
    <td><a href=\"modelica://OpenTEMPEST.SOC.Stack\">SOC.Stack</a></td></tr>
<tr><td>4. Module</td>
    <td>Several stacks in a shared thermal environment (examples in preparation).</td>
    <td><a href=\"modelica://OpenTEMPEST.SOC.MultiReactor\">SOC.MultiReactor</a></td></tr>
<tr><td>5. System</td>
    <td>Stacks combined with balance of plant: heat exchangers, heaters, blowers, compressors, reformers,
    afterburners, downstream reactors, power electronics and controllers.</td>
    <td><a href=\"modelica://OpenTEMPEST.BOP\">BOP</a>,
        <a href=\"modelica://OpenTEMPEST.Flow\">Flow</a>,
        <a href=\"modelica://OpenTEMPEST.Blocks\">Blocks</a></td></tr>
</table>
<h4>How the levels are connected</h4>
<p>A stack contains an array of cells <code>cell[Ncell]</code>, a cell contains a PEN, and the PEN contains one
electrochemistry object per control volume <code>electrochem[N]</code>. The electrochemistry is a
<code>replaceable model Electrochem</code>: you choose it once at the level you work on, and every level passes it
down with <code>redeclare model Electrochem = Electrochem</code>. Gas media (<code>FCMedium</code>,
<code>ACMedium</code>) and the flux interpolation scheme (<code>fluxInterp</code>) are passed down the same way.</p>
<p>Supporting packages:</p>
<ul>
<li><a href=\"modelica://OpenTEMPEST.Heat\">Heat</a>: 0D to 3D solids, conduction, convection, radiation and heat port adaptors.</li>
<li><a href=\"modelica://OpenTEMPEST.Medium\">Medium</a>: ideal gas mixtures for fuel, air and exhaust.</li>
<li><a href=\"modelica://OpenTEMPEST.Solid\">Solid</a>: temperature dependent properties of solid materials.</li>
<li><a href=\"modelica://OpenTEMPEST.Units\">Units</a> and <a href=\"modelica://OpenTEMPEST.Enumerations\">Enumerations</a>: types used across the library.</li>
</ul>
</html>"));
    end Overview;

    class GettingStarted "Getting started: run a first simulation"
      extends Modelica.Icons.Information;
      annotation (Documentation(info="<html>
<h4>Requirements</h4>
<ul>
<li>Dymola (tested). OpenModelica works for many models but is not tested systematically yet.</li>
<li>Modelica Standard Library 3.2.3.</li>
<li><a href=\"https://github.com/casella/ThermoPower\">ThermoPower</a> 3.1.</li>
</ul>
<p>The migration to MSL 4.x without ThermoPower is in progress. See the
<a href=\"modelica://OpenTEMPEST.UsersGuide.ReleaseNotes\">release notes</a>.</p>
<h4>First simulation</h4>
<ol>
<li>Load the Modelica Standard Library, ThermoPower and OpenTEMPEST (open their <code>package.mo</code> files).</li>
<li>Open <a href=\"modelica://OpenTEMPEST.Examples.Stack.Stack1DTest\">Examples.Stack.Stack1DTest</a>.
    It connects a stack of 1D co-flow cells to gas sources and sinks and to a current controlled power supply,
    and runs steam electrolysis with a stepwise current ramp.</li>
<li>Simulate it. The stop time and the solver are stored in the model.</li>
<li>Plot for example <code>stack.cell[3].pen.T</code> (PEN temperature along the cell),
    <code>stack.cell[3].pen.J</code> (current density) and <code>stack.cell[3].pen.Uop[1]</code> (cell voltage).</li>
</ol>
<p>Faster alternatives are <a href=\"modelica://OpenTEMPEST.Examples.Stack.Stack0DTest\">Stack0DTest</a>
(lumped stack) and <a href=\"modelica://OpenTEMPEST.Examples.Stack.SimplifiedStackTest\">SimplifiedStackTest</a>
(a few detailed cells, the rest simplified).</p>
<h4>Next steps</h4>
<p>The user workflows in the repository (<code>docs/UserWorkflows.md</code>) explain step by step how to build a
cell test bench, choose or write an electrochemistry model, work with stack models, build a system model and
extend the library in your own package. The Markdown library reference (<code>docs/library-reference</code>)
lists every class with its parameters and connectors.</p>
<h4>Practical tips</h4>
<ul>
<li>Keep your own models in your own package and use OpenTEMPEST models through modifiers, <code>extends</code>
    and <code>redeclare</code>, instead of editing the library.</li>
<li>Start from zero current (<code>Jstart=0</code>) and ramp the current up once the temperatures have settled.</li>
<li>Set the start values in the <em>Initialization</em> tab close to the expected operating point.
    Poor start values are the most common reason for a failed initialisation.</li>
<li>Use a stiff solver such as Cvode or Dassl.</li>
</ul>
</html>"));
    end GettingStarted;

    class Conventions "Conventions used in the library"
      extends Modelica.Icons.Information;
      annotation (Documentation(info="<html>
<h4>Signs</h4>
<ul>
<li>The current density <code>J</code> is positive in fuel cell (SOFC) mode and negative in electrolysis
    (SOEC) mode.</li>
<li>Heat flows into a component through its heat ports are positive.</li>
<li>Mass flow rates at a connector are positive when the gas flows into the component.</li>
</ul>
<h4>Gas composition</h4>
<p>The order of species follows the medium. The index into mass fractions <code>X</code> and mole fractions
<code>y</code> depends on it:</p>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"4\">
<tr><th>Medium</th><th>Species, in index order</th></tr>
<tr><td><a href=\"modelica://OpenTEMPEST.Medium.Fuel_CH4\">Fuel_CH4</a></td><td>H2, CH4, CO2, CO, H2O, N2</td></tr>
<tr><td><a href=\"modelica://OpenTEMPEST.Medium.Fuel_C3\">Fuel_C3</a></td><td>H2, CH4, CO2, CO, H2O, N2, C3H8</td></tr>
<tr><td><a href=\"modelica://OpenTEMPEST.Medium.Fuel_MethanolReactor\">Fuel_MethanolReactor</a></td><td>H2, CH4, CO2, CO, H2O, N2, CH3OH</td></tr>
<tr><td><a href=\"modelica://OpenTEMPEST.Medium.Air_Medium\">Air_Medium</a></td><td>O2, N2</td></tr>
<tr><td><a href=\"modelica://OpenTEMPEST.Medium.SimpleExhaust\">SimpleExhaust</a></td><td>O2, CO2, H2O, N2</td></tr>
<tr><td><a href=\"modelica://OpenTEMPEST.Medium.FlueGas\">FlueGas</a></td><td>O2, Ar, H2O, CO2, N2</td></tr>
</table>
<h4>Geometry and discretisation</h4>
<ul>
<li><code>x</code> is the fuel flow direction, <code>y</code> the cell width (the air flow direction in cross-flow
    cells) and <code>z</code> the stacking direction.</li>
<li><code>lX</code>, <code>lY</code>, <code>lZ</code> are lengths in these directions; <code>N</code>, <code>nX</code>,
    <code>nY</code> are numbers of control volumes; <code>Ncell</code> is the number of cells in a stack.</li>
<li>Abbreviations: PEN (positive electrode, electrolyte, negative electrode), FC (fuel channel), AC (air channel),
    IC (interconnect), ASR (area specific resistance), BoP (balance of plant).</li>
</ul>
<h4>Units</h4>
<p>All variables are in SI units. Correlations fitted in other units say so in their description, for example the
ASR correlations in <a href=\"modelica://OpenTEMPEST.SOC.Electrochem.ASR\">SOC.Electrochem.ASR</a>, which are
written for Ohm.cm<sup>2</sup> but return Ohm.m<sup>2</sup>.</p>
<h4>Documentation</h4>
<p>Every class and parameter has a short description string (sentence case, no full stop). Longer documentation is
in the info text, with sections as <code>&lt;h4&gt;</code> headings. The rules are written down in
<code>docs/DocumentationStyle.md</code> in the repository.</p>
</html>"));
    end Conventions;

    class ReleaseNotes "Release notes"
      extends Modelica.Icons.ReleaseNotes;
      annotation (Documentation(info="<html>
<h4>Development version</h4>
<ul>
<li>Uniform description strings for all classes and info texts following one style.</li>
<li>User's guide in the library, user workflows and a Markdown library reference in the repository.</li>
<li>Planned: migration to the Modelica Standard Library 4.x and own gas and heat interfaces, so that
    ThermoPower is no longer needed.</li>
</ul>
<h4>Version 0.1.0</h4>
<p>First public release of OpenTEMPEST, based on the Modelica Standard Library 3.2.3 and ThermoPower 3.1.</p>
<ul>
<li>Electrochemistry: lumped ASR and detailed Butler-Volmer models for steam, CO2 and co-electrolysis.</li>
<li>Cells: 1D co-flow and 2D cross-flow cells, and a cross-flow cell with simplified 1D channels.</li>
<li>Stacks: 0D, 1D and cross-flow stacks, and simplified stacks with a few detailed cells.</li>
<li>Balance of plant: heat exchangers, heaters, blowers, compressors, reformers, afterburner,
    methanol synthesis and power electronics.</li>
<li>Examples for the stack models.</li>
</ul>
</html>"));
    end ReleaseNotes;

    class References "Publications"
      extends Modelica.Icons.References;
      annotation (Documentation(info="<html>
<p>If you use OpenTEMPEST in your work, please cite the software
(<a href=\"https://doi.org/10.5281/zenodo.21771899\">doi:10.5281/zenodo.21771899</a>) and the publication that
describes the models you used.</p>
<ul>
<li>J. Hollmann et al. 2026. Improving low-load efficiency of marine SOFC powertrains by means of efficient
    module hot standby operation. <em>J. Electrochem. Soc.</em> 173, 124503.
    <a href=\"https://doi.org/10.1149/1945-7111/ae7479\">doi:10.1149/1945-7111/ae7479</a></li>
<li>R. Lorenz et al. 2025. Control concept for solid oxide electrolysis reactors to promote highly transient
    operation in modular plants. <em>Energy Convers. Manag. X</em> 28, 101236.
    <a href=\"https://doi.org/10.1016/j.ecmx.2025.101236\">doi:10.1016/j.ecmx.2025.101236</a></li>
<li>F. Sedeqi et al. 2024. Solid oxide cell reactor model for transient and stationary electrochemical H2O and
    CO2 conversion process studies. <em>J. Electrochem. Soc.</em> 171, 074507.
    <a href=\"https://doi.org/10.1149/1945-7111/ad5e01\">doi:10.1149/1945-7111/ad5e01</a></li>
<li>M. Tomberg et al. 2023. Operation strategies for a flexible megawatt scale electrolysis system for synthesis
    gas and hydrogen production with direct air capture of carbon dioxide. <em>Sustainable Energy Fuels</em> 7(2),
    471-484. <a href=\"https://doi.org/10.1039/d2se01473d\">doi:10.1039/d2se01473d</a></li>
<li>M. Tomberg et al. 2022. A new approach to modeling solid oxide cell reactors with multiple stacks for process
    system simulation. <em>J. Electrochem. Soc.</em> 169, 054530.
    <a href=\"https://doi.org/10.1149/1945-7111/ac7009\">doi:10.1149/1945-7111/ac7009</a></li>
<li>M. Tomberg et al. 2022. Operation of a solid oxide fuel cell reactor with multiple stacks in a pressured system
    with fuel gas recirculation. <em>Energy Technol.</em> 10.
    <a href=\"https://doi.org/10.1002/ente.202101075\">doi:10.1002/ente.202101075</a></li>
<li>S. Srikanth et al. 2018. Transient reversible solid oxide cell reactor operation: experimentally validated
    modeling and analysis. <em>Appl. Energy</em> 232, 473-488. <a href=\"https://doi.org/10.1016/j.apenergy.2018.09.186\">doi:10.1016/j.apenergy.2018.09.186</a></li>
</ul>
</html>"));
    end References;

    class Contact "Contact and license"
      extends Modelica.Icons.Contact;
      annotation (Documentation(info="<html>
<h4>Contact</h4>
<p>OpenTEMPEST is developed at the German Aerospace Center (DLR), Institute of Engineering Thermodynamics.</p>
<ul>
<li>E-mail: <a href=\"mailto:tempest@dlr.de\">tempest@dlr.de</a></li>
<li>Issues and merge requests:
    <a href=\"https://codeberg.org/DLR-TT-ESI-EHT/OpenTEMPEST\">codeberg.org/DLR-TT-ESI-EHT/OpenTEMPEST</a></li>
</ul>
<p>Feedback, bug reports and contributions are welcome. Please follow the documentation conventions in
<a href=\"modelica://OpenTEMPEST.UsersGuide.Conventions\">Conventions</a> when you contribute a model.</p>
<h4>License</h4>
<p>Copyright &copy; 2025, DLR-TT. OpenTEMPEST is released under the BSD 3-Clause License. The full text is in the
<code>LICENSE</code> file of the repository.</p>
</html>"));
    end Contact;

    annotation (DocumentationClass=true, Documentation(info="<html>
<p>This guide explains how OpenTEMPEST is organised, how to run a first simulation and which conventions the
library follows.</p>
<ul>
<li><a href=\"modelica://OpenTEMPEST.UsersGuide.Overview\">Overview</a> of the modelling levels</li>
<li><a href=\"modelica://OpenTEMPEST.UsersGuide.GettingStarted\">Getting started</a></li>
<li><a href=\"modelica://OpenTEMPEST.UsersGuide.Conventions\">Conventions</a></li>
<li><a href=\"modelica://OpenTEMPEST.UsersGuide.ReleaseNotes\">Release notes</a></li>
<li><a href=\"modelica://OpenTEMPEST.UsersGuide.References\">Publications</a></li>
<li><a href=\"modelica://OpenTEMPEST.UsersGuide.Contact\">Contact and license</a></li>
</ul>
</html>"));
  end UsersGuide;

  annotation (
    preferredView="info",
    uses(
      Modelica(version="3.2.3"),
      ThermoPower(version="3.1"),
      ModelicaServices(version="3.2.3")),
    version="0.1.0",
    Documentation(info="<html>
<p><strong>OpenTEMPEST</strong> (<strong>T</strong>ransient <strong>E</strong>lectrochemical <strong>M</strong>odel
for <strong>P</strong>rocess and <strong>E</strong>nergy <strong>S</strong>ys<strong>T</strong>ems) is an open source
Modelica library for the dynamic simulation of electrochemical energy systems, from a single cell to a complete plant.
It covers high temperature solid oxide cells in fuel cell and electrolysis mode (SOFC/SOEC), and alkaline
electrolysis (AEL) is in preparation.</p>
<p>The library is built for steady-state and transient studies at process system level: part load and hot standby
operation, start-up and shutdown, control concepts for fluctuating power, and the coupling of electrolysers or fuel
cells to upstream and downstream processes. Validated simplifications at every modelling level let you scale from
the cell to stacks and multi-stack modules without losing the thermal behaviour that limits real operation.</p>
<h4>Highlights</h4>
<ul>
<li>Finite volume cell models in 1D (co-flow) and 2D (cross-flow) with PEN, gas channels and interconnects,
    including internal reforming and water gas shift in the fuel channel.</li>
<li>Replaceable electrochemistry: lumped ASR models or detailed Butler-Volmer models with ohmic and diffusion
    losses for steam, CO2 and co-electrolysis.</li>
<li>Stack models from a lumped 0D stack to stacks of 2D cells, and simplified stacks that resolve a few cells in
    detail and lump the rest.</li>
<li>Balance of plant components with off-design behaviour: heat exchangers, electric heaters, blowers,
    compressors, pre-reformers, afterburners, methanol synthesis and power electronics.</li>
</ul>
<h4>Package structure</h4>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"4\">
<tr><th>Package</th><th>Contents</th></tr>
<tr><td><a href=\"modelica://OpenTEMPEST.UsersGuide\">UsersGuide</a></td><td>Overview, getting started, conventions, release notes, publications, contact and license</td></tr>
<tr><td><a href=\"modelica://OpenTEMPEST.Examples\">Examples</a></td><td>Ready to run example models</td></tr>
<tr><td><a href=\"modelica://OpenTEMPEST.SOC\">SOC</a></td><td>Solid oxide cell models: electrochemistry, cells and stacks</td></tr>
<tr><td><a href=\"modelica://OpenTEMPEST.BOP\">BOP</a></td><td>Balance of plant components</td></tr>
<tr><td><a href=\"modelica://OpenTEMPEST.Flow\">Flow</a></td><td>Gas flow components: sources, splitters, manifolds and sensors</td></tr>
<tr><td><a href=\"modelica://OpenTEMPEST.Heat\">Heat</a></td><td>Solids and heat transfer</td></tr>
<tr><td><a href=\"modelica://OpenTEMPEST.Blocks\">Blocks</a></td><td>Signal sources, controllers and helper functions</td></tr>
<tr><td><a href=\"modelica://OpenTEMPEST.Medium\">Medium</a></td><td>Gas mixture media</td></tr>
<tr><td><a href=\"modelica://OpenTEMPEST.Solid\">Solid</a></td><td>Solid material properties</td></tr>
<tr><td><a href=\"modelica://OpenTEMPEST.AEL\">AEL</a></td><td>Alkaline electrolysis (in preparation)</td></tr>
<tr><td><a href=\"modelica://OpenTEMPEST.Units\">Units</a>, <a href=\"modelica://OpenTEMPEST.Enumerations\">Enumerations</a></td><td>Types used across the library</td></tr>
</table>
<h4>Where to start</h4>
<p>Read the <a href=\"modelica://OpenTEMPEST.UsersGuide.Overview\">overview</a> of the modelling levels, then
follow <a href=\"modelica://OpenTEMPEST.UsersGuide.GettingStarted\">getting started</a> to run
<a href=\"modelica://OpenTEMPEST.Examples.Stack.Stack1DTest\">Stack1DTest</a>.</p>
<p>Copyright &copy; 2025, DLR-TT. Released under the BSD 3-Clause License, see
<a href=\"modelica://OpenTEMPEST.UsersGuide.Contact\">Contact and license</a>.</p>
</html>"));
end OpenTEMPEST;
