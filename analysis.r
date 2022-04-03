if (!require("R.matlab", quietly = TRUE)){install.packages("R.matlab")}
if (!require("tidyverse", quietly = TRUE)){install.packages("tidyverse")}

library(R.matlab)
library(tidyverse)

ecgsamp <- readMat("ecgsamp.mat")
ecg <- ecgsamp$ecg
rm(ecgsamp)

mu <- mean(ecg)
sigma <- sd(ecg)

ggplot() + geom_line(aes(x = 1:length(ecg)/250, y = ecg*0.001)) + geom_hline(yintercept = mu*0.001 + 2.576*sigma*0.001, color = 'red', linetype = 'dashed') + xlim(120, 180) + labs(x = 'Time (s)', y = 'Potential (mV)', title = 'ECG from 120s to 180s, Decision Line') + theme_classic()

peaking_times <- which(ecg > mu + 2.576*sigma)

peak_times <- numeric(0)
prev_time <- peaking_times[1]
prev_peak <- peaking_times[1]
for (time in peaking_times[-1]) {
  if (time == prev_time + 1) {
    if (ecg[time] > ecg[prev_peak]) {
      prev_peak <- time
    }
    if (time == peaking_times[length(peaking_times)]) {
      peak_times <- append(peak_times, prev_peak)
    }
  } else {
    peak_times <- append(peak_times, prev_peak)
    prev_peak <- time
  }
  prev_time <- time
}

ggplot() + geom_line(aes(x = 1:length(ecg)/250, y = ecg*0.001)) + geom_point(aes(x = peak_times/250, y = ecg[peak_times]*0.001), color = 'red') + xlim(120, 180) + geom_hline(yintercept = mu*0.001 + 2.576*sigma*0.001, color = 'red', linetype = 'dashed') + labs(title = 'ECG from 120s to 180s (Highlighted Peaks, Decision Line)', x = 'Time (s)', y = 'Potential (mV)') + theme_classic()



rmssd_series <- function(time.series, step, epoch_width) {
  rmssd.series <- numeric((time.series[length(time.series)] - epoch_width) %/% step)
  rmssd.time <- numeric((time.series[length(time.series)] - epoch_width) %/% step)
  rmssd.n <- numeric((time.series[length(time.series)] - epoch_width) %/% step)
  for (i in 0:(((time.series[length(time.series)] - epoch_width) %/% step) - 1)) {
    times <- time.series[time.series >= step*i & time.series <= (step*i + epoch_width)]
    n <- length(times) - 1
    if (n < 2) {
      return(c('Error, epoch width too small.', times))
    }
    ssds <- numeric(n-1)
    for (j in 1:(n-1)) {
      ssds[j] <- ((times[j+2] - times[j+1]) - (times[j+1] - times[j]))^2
    }
    ssd <- sum(ssds)
    rmssd <- sqrt(1/(n-1)*ssd)
    rmssd.series[i+1] <- rmssd
    rmssd.time[i+1] <- step*i + epoch_width/2
    rmssd.n <- n
  }
  rmssd.frame <- data.frame(time = rmssd.time, rmssd = rmssd.series, n = rmssd.n)
  return(rmssd.frame)
}

rmssd.data <- rmssd_series(peak_times, step = 500, epoch_width = 250*30)
rmssd.data %>% ggplot(aes(x = time/250, y = log(rmssd))) + geom_line() + labs(title = 'RMSSD for ECG (30s intervals)', x = 'Time (s)', y = 'Log of RMSSD') + scale_x_continuous(breaks = seq(30, 360, by = 30)) + theme_classic()

rmssd.data <- rmssd_series(peak_times, step = 500, epoch_width = 250*240)
rmssd.data %>% ggplot(aes(x = time/250, y = log(rmssd))) + geom_line() + labs(title = 'RMSSD for ECG (240s intervals)', x = 'Time (s)', y = 'Log of RMSSD') + scale_x_continuous(breaks = seq(30, 360, by = 30)) + theme_classic()
