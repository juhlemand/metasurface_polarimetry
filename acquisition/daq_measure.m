function data = daq_measure(duration, filename)
% This function compels the daq unit to take a measurement.
  global adc

  adc.DurationInSeconds = duration;
  data = startForeground(adc);
  csvwrite(filename, data);
  data = mean(data);
  adc.release();
end

