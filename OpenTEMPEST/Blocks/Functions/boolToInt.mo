within OpenTEMPEST.Blocks.Functions;
function boolToInt "return 1 for true and 0 for false entries"
  extends Modelica.Icons.Function;

  input Boolean boolVec[:] "Boolean vector";
  output Integer intVec[size(boolVec, 1)]
      "integer vector";

algorithm
    for i in 1:size(boolVec, 1) loop
      if boolVec[i] then
        intVec[i] := 1;
      else
        intVec[i] := 0;
      end if;
    end for;

end boolToInt;
