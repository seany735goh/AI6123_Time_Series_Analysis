y = read.delim("Project_1_Data.txt", header = TRUE)
y = y$x
x = length(y)
plot(1:x,y)
lines(1:x,y, type="l")
acf(y, lag.max=50)
pacf(y, lag.max=50)

# Differencing once, plot time series, acf and pacf
yd1 = diff(y)
plot(1:(x-1),yd1)
lines(1:(x-1),yd1, type="l")
acf(yd1, lag.max=50)
pacf(yd1, lag.max=50)

# Fit MA model
fit1ma = arima(y, order=c(0,1,24))
tsdiag(fit1ma)
fit1ma

# Fit AR model
fit1ar = arima(y, order=c(3,1,0))
tsdiag(fit1ar)
fit1ar

# Use auto.ARIMA
fit1auto = auto.arima(y, d=1, max.p=30, max.q=30, ic="aic")
tsdiag(fit1auto) # Produces ARIMA(1,1,1)
fit1auto

# Differencing twice, plot time series, acf and pacf
yd2 = diff(yd1)
plot(1:(x-2),yd2)
lines(1:(x-2),yd2, type="l")
acf(yd2, lag.max=50)
pacf(yd2, lag.max=50)

# Fit MA model
fit2ma = arima(y, order=c(0,2,3))
tsdiag(fit2ma)
fit2ma

# Fit AR model
fit2ar = arima(y, order=c(2,2,0))
tsdiag(fit2ar)
fit2ar

# Use auto.ARIMA
fit2auto = auto.arima(y, d=2, max.p=30, max.q=30, ic="aic")
tsdiag(fit2auto) # Produces ARIMA(2,2,0)
fit2auto

# Plot Time-Series Predictions of ARIMA Models
plot(1:x, y, main="Time-Series Predictions of ARIMA Models")
lines(1:x, y, type="l")
lines(1:x, y-fit1ma$residuals, type="l", col="red")
lines(1:x, y-fit1ar$residuals, type="l", col="blue")
lines(1:x, y-fit1auto$residuals, type="l", col="purple")
lines(1:x, y-fit2ma$residuals, type="l", col="green")
lines(1:x, y-fit2ar$residuals, type="l", col="orange")
legend(1, 225, legend=c("ARIMA(0,1,24)", "ARIMA(3,1,0)", "ARIMA(1,1,1)", "ARIMA(0,2,3)", "ARIMA(2,2,0)"), col=c("red", "blue", "purple", "green", "orange"), lty=1:1)
