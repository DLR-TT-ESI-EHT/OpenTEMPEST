within OpenTEMPEST.Blocks.Functions;
function limit
  "Limits the input variable between two boundaries"
  extends Modelica.Icons.Function;
  input Real minVal;
  input Real u;
  input Real maxVal;
  output Real y;
algorithm
  y := min(maxVal,max(minVal,u));

  annotation (
    Documentation(info="<html>
This function limits the input variable between two boundaries. 
</html>",
        revisions="<html>
</html>"));
end limit;
