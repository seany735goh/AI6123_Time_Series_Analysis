# This project is to analyze financial data. 
# The data are from the daily historical Apple stock 
# prices(open, high, low, close and adjusted prices) 
# from February 1, 2002 to January 31, 2017 extracted from 
# the Yahoo Finance website. The data has logged the prices of the 
# Apple stock everyday and comprises of the open, close, low, high and the 
# adjusted close prices of the stock for the span of 15 years.
# The goal of the project is to discover an interesting trend in the apple 
# stock prices over the past 15 years (3775 attributes) and to design and 
# develop the best model for forecasting.
# Source:    http://finance.yahoo.com/quote/AAPL/history?p=AAPL

##### Install/import libraries

install.packages("forecast") # for BoxCox function
install.packages("lubridate") # for HW functions
install.packages("TSA")
install.packages("timeDate")
install.packages("changepoint")
install.packages("fGarch")
install.packages("rugarch")
install.packages("zoo")
library(forecast)
library(lubridate)
library(tseries)
library(TSA)
library(timeDate)
library(changepoint)
library(fGarch)
library(rugarch)
library(zoo)
library(quantmod) #getSymbols

##### Import data

y = getSymbols('AAPL', from='2002-02-01', to='2017-02-01', src='yahoo', auto.assign = F) 
y = na.omit(y)
y_close = y[,4]

# y2 = getSymbols('AAPL', from='2017-02-01', to='2017-04-01', src='yahoo', auto.assign = F) 
y2 = getSymbols('AAPL', from='2017-02-01', to='2018-02-01', src='yahoo', auto.assign = F) 
y2 = na.omit(y2)
y_close2 = y2[,4]

y_close_c = c(y_close, y_close2)

dates = index(y) # 3776 datapoints
dates2 = index(y2) # 252 datapoints
dates_c = c(dates, dates2)

plot(dates, y_close, ylim=c(0,35), xlab='Year', ylab='USD', main='Apple Stock Close Price', type='l')
par(mfrow=c(2,1))
acf(y_close, main='Close Price')
pacf(y_close, main='Close Price')
adf.test(y_close) #0.4114, therefore reject stationary hypothesis, therefore is non-stationary

##### Preprocess ARIMA

lambda = BoxCox.lambda(y_close) # 0.114639253068155
y_transformed = BoxCox(y_close, lambda)

par(mfrow=c(1,1))
plot(dates, y_transformed, ylim=c(-2,5), xlab='Year', ylab='Transformed USD', main='Transformed Close Price', type='l')
par(mfrow=c(2,1))
acf(y_transformed, main='Transformed Close Price')
pacf(y_transformed, main='Transformed Close Price')

y_d = diff(y_transformed)
y_d = na.omit(y_d)
par(mfrow=c(1,1))
plot(dates[1:3775], y_d, xlab='Year', ylab='Differenced Transformed USD', main='Differenced Transformed Close Price', type='l')
adf.test(y_d)
kpss.test(y_d)

par(mfrow=c(2,1))
acf(y_d, main='Differenced Transformed Close Price') # q=4
pacf(y_d, main='Differenced Transformed Close Price') # p=4

##### ARIMA fitting

fit_ar = arima(y_transformed, order=c(0,1,4))
# fit_ar = auto.arima(y_transformed, d=1, max.p=30, max.q=30, ic="aic") # (2,1,2)
tsdiag(fit_ar)
acf(fit_ar$residuals, lag.max=50)
pacf(fit_ar$residuals, lag.max=50)
AIC(fit_ar) 
# -16839.62 for (4,1,0)
# -16839.49 for (0,1,4)
# -16827.33 for (16,1,4)
# -16842.67 for (4,1,16)
# -16845.84 for (4,1,4)

par(mfrow=c(1,1))
plot(1:length(dates), y_transformed, xlim=c(0,4200), ylim=c(-2,5), type='l')
forecast_fit_ar = predict(fit_ar, n.ahead=length(dates2))
lines(((length(dates)+1):(length(dates)+length(dates2))), forecast_fit_ar$pred, col='red')

forecast_y = InvBoxCox(forecast_fit_ar$pred, lambda)
forecast_ub = InvBoxCox(forecast_fit_ar$pred + 1.96*forecast_fit_ar$se, lambda) # upper bound
forecast_lb = InvBoxCox(forecast_fit_ar$pred - 1.96*forecast_fit_ar$se, lambda) # lower bound

plot(dates_c, y_close_c, ylim=c(0,50), xlab='Year', ylab='USD', main='Apple Stock Close Price', type='l')
lines(dates2, forecast_ub, type='l', col='blue')
lines(dates2, forecast_lb, type='l', col='blue')
lines(dates2, forecast_y, type='l', col='red')
lines(c(dates_c[3777],dates_c[3777]), c(0,50), lty=2, col='black')
legend("topleft", 
       c("Close Price Data (01/02/2002 to 31/01/2018)", 
         "Model Forecast (01/02/2017 to 31/01/2018)",  
         "Forecast Upper & Lower Bounds"), 
       col=c('black', 'blue', 'red'), 
       lty=c(1,1,1,1))

##### GARCH Analysis
# returns = diff(y_transformed)*100
# returns = diff(log(y_close))*100 # log transformation, differencing, then multiplied by 100
lambda = 0
returns = diff(BoxCox(y_close, lambda))*100 # log transformation, differencing, then multiplied by 100
returns_2 = diff(BoxCox(y_close2, lambda))*100 # log transformation, differencing, then multiplied by 100
returns_c = diff(BoxCox(c(y_close, y_close2), lambda))*100 # log transformation, differencing, then multiplied by 100
returns = returns[!is.na(returns)] 
par(mfrow=c(1,1))
plot(dates[1:3775], returns, xlab='Year', ylab='Percentage Returns', main='Percentage Returns', type='l')
abline(h=0, col='red')

par(mfrow=c(2,1))
acf(returns) # q = 2
pacf(returns) # p = 2
adf.test(returns) # p-value = 0.01, therefore reject h0 (non-stationary), therefore is stationary
kpss.test(returns) # p-value = 0.1, means fail to reject h0 (level stationary), therefore is level stationary

acf(abs(returns))
pacf(abs(returns))
# Significant autocorrelations observed, hence returns are not i.i.d

acf(returns^2)
pacf(returns^2)
# Significant autocorrelations observed, hence returns are not i.i.d


par(mfrow=c(1,1))
qqnorm(returns, main='Normal Q-Q Plot of Returns'); qqline(returns) 
# data varies largely away from qqline, therefore heavy tailed distribution
kurtosis(returns) # 4.915312 (BoxCox) OR 5.440092 (log), therefore heavy tailed distribution
# ?skewness

##### EACF Observation

eacf(returns) #GARCH(0,2) 4x in triangle, OR GARCH(0,4) 0x in triangle
eacf(abs(returns)) #GARCH(1,1) 5x in triangle
eacf(returns^2) #GARCH(2,2) 14x in triangle, OR GARCH(1,3) 15x in triangle

##### GARCH FITTING

g1 = garch(returns, order=c(0,2)) 
summary(g1) # Box-Ljung p-value = 0.4253
AIC(g1) # 16547.11

g2 = garch(returns, order=c(0,4))
summary(g2) # Box-Ljung p-value = 0.3304
AIC(g2) # 16447.81

g3 = garch(returns, order=c(1,1))
summary(g3) # Box-Ljung p-value = 0.2878
AIC(g3) # 16205.51

g4 = garch(returns, order=c(1,8))
summary(g4) # Box-Ljung p-value = 0.8546
AIC(g4) # 16313.33

g5 = garch(returns, order=c(1,10))
summary(g5) # Box-Ljung p-value = 0.838
AIC(g5) # 16272.13

##### Fitted GARCH Check

plot(residuals(g3), type='h', ylab='Standardized Residuals')
qqnorm(residuals(g3),  main='Normal Q-Q Plot of Standardized Residuals from Fitted GARCH(1,1) Model'); qqline(residuals(g3))
kurtosis(residuals(g3)) # 3.540995 for log

acf(residuals(g3)^2, na.action=na.omit) # squared residuals are serially uncorrelated
gBox(g3, method='squared')
logLik(g3) # -8099.754

##### Forecasting Fit

s = ugarchspec(mean.model = list(armaOrder = c(0,0)),
               variance.model = list(model = 'eGARCH', garchOrder=c(1,1)), #sGARCH, gjrGARCH, eGARCH
               distribution.model = 'std') # norm, snorm, std, sstd
fit_garch = ugarchfit(data=returns, spec=s)
fit_garch
# Model, distribution, loglik, AIC
# sGARCH, norm, -8088.507, 4.2874
# sGARCH, snorm, -8088.274, 4.2878
# sGARCH, std, -7912.203, 4.1945
# sGARCH, sstd, -7911.74, 4.1948
# gjrGARCH, std, -7902.228, 4.1898
# eGARCH, std, -7888.421, 4.1825

plot(fit_garch, which=9) #9, 10, 11

setfixed(s) = as.list(coef(fit_garch))
forecast_garch = ugarchforecast(data=returns, fitORspec = s, n.ahead = 365) #n.ahead includes weekends, therefore 365

# plot forecasted returns
plot(forecast_garch, which=1)
lines(dates_c[3776:4028], returns_c[3776:4028], col='black')
legend("topleft", 
       c("Returns Data (07/09/2016 - 31/01/2017)", 
         "Forecast (01/02/2017 - 31/01/2018)",
         "Returns Data (01/02/2017 - 31/01/2018)"),
       col=c('blue', 'red', 'black'), 
       lty=c(1,1,1))

# plot forecasted variance
plot(forecast_garch, which=3)
plot(sigma(forecast_garch)/100)

# plot Forcasted close price
count=0
max_value = 0
while (max_value < 36){
  forecast_returns = ugarchpath(spec = s, m.sim = 3, n.sim = 365)
  forecast_close = exp(apply(fitted(forecast_returns)/100, 2, "cumsum")) + as.numeric(y_close[3776])
  count = count + 1
  max_value = max(forecast_close)
  if (count%%10 == 0){
    print(count/10)
  } 
}
matplot(forecast_close, type='l')
  
dates_forecast_close = seq(as.Date("2017/02/01"), by = "day", length.out = 365)
plot(dates_c, y_close_c, ylim=c(0,50), xlab='Year', ylab='USD', main='Apple Stock Close Price', type='l')
lines(dates_forecast_close, forecast_close[,1], col='red')
lines(dates_forecast_close, forecast_close[,2], col='blue')
lines(dates_forecast_close, forecast_close[,3], col='orange')
lines(c(dates_c[3777],dates_c[3777]), c(0,50), lty=2, col='black')
legend("topleft", 
       c("Close Price Data (01/02/2002 - 31/01/2018)", 
         "Simulated Forecast 1 (01/02/2017 - 31/01/2018)", 
         "Simulated Forecast 2 (01/02/2017 - 31/01/2018)", 
         "Simulated Forecast 3 (01/02/2017 - 31/01/2018)" ),
       col=c('black', 'red', 'blue', 'orange'), 
       lty=c(1,1,1,1))

# Plot zoomed in close price on 01/02/2017

plot(dates_c[3677:3877], y_close_c[3677:3877], ylim=c(25,40), xlab='Month', ylab='USD', main='Apple Stock Close Price', type='l')
lines(c(dates_c[3777],dates_c[3777]), c(25,40), lty=2, col='red')
legend("topleft", 
       c("Close Price Data (08/09/2016 - 26/06/2017)", 
         "01/02/2017"),
       col=c('black', 'red'), 
       lty=c(1,1))

# Detect change points

cpt_pelt = cpt.mean(y_close_c_num, penalty='MBIC', method='PELT', minseglen=5)
plot(cpt_pelt, type='l', cpt.col='red', cpt.width=3, xlab='Index', ylab='USD', main='Apple Stock Close Price')
lines(c(3777,3777), c(0,50), lty=2)
moo = data.frame(cpt_pelt@cpts, cpt_pelt@param.est[['mean']])

