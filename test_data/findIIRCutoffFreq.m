
%finds approximate cutoff frequency using eq 8-12 pg 90 of DSP filters book
%first find f1 and then uses the relationship with Q and f0 to find f2
function [f1, f2] = findIIRCutoffFreq(fs, f0, Q)

  theta0 = 2 * pi * f0 / fs;
  beginTheta = 0;
  endTheta = theta0;
  if (theta0 > 0.01)
    step = 0.01;
  elseif (theta0 > 1e-4)
    step = 1e-4;
  else
    step = 1e-6;
  endif

  while(step > 1e-16)
    x = beginTheta : step : endTheta;
    result = (sin(x) .* tan(theta0 ./ (2 .* Q))) ./ sqrt((sin(x) .* tan(theta0 ./ (2 .* Q))).^2 + (cos(x) - cos(theta0)).^2) - 1/sqrt(2);
    %indicate zero crossing with sign change
    signResult = sign(result);

    %technically possible we could have a step land right on zero so check for it
    if (isempty(find(signResult == 0)) == false)
      beginTheta = (find(signResult == 0, 1) - 1) * step + beginTheta;
      break;
    endif

    endTheta = (find(signResult == 1, 1) - 1) * step + beginTheta;
    beginTheta = endTheta - step;
    step = step * 0.01;
  endwhile

  f1 = beginTheta * fs / (2 * pi);
  f2 = f0 / Q + f1;

endfunction
