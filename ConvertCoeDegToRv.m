function [rVecECI, vVecECI] = ConvertCoeDegToRv(coe)

[rVecECI, vVecECI] = ConvertCoeToRv([coe(1:2) deg2rad(coe(3:6))]);

end
