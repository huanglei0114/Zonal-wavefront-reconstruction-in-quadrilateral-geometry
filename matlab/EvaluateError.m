function [Error_wo_mean, rms_Error_wo_mean, Error] = EvaluateError(Eval, True)

Error = Eval - True;
Error_wo_mean = Error - mean(Error(isfinite(Error)));
rms_Error_wo_mean = rms(Error_wo_mean(isfinite(Error_wo_mean)));

end
