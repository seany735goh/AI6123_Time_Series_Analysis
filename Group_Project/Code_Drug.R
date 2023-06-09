# # Packages
install.packages("forecast") # for BoxCox function
install.packages("lubridate") # for HW functions
library(forecast)
library(lubridate)
library(tseries)
# setwd('C:/Users/guang/OneDrive/Desktop/AI6123/Project 2/repo')

y = read.delim("drug.txt", sep = ",", header = TRUE)
# Monthly data --> s = 12, increasing trend
y = y$value
plot(y, main='Time Series: Drug') 
lines(y, type="l") 
# Increasing variance over time, thus use boxcox transformation

lambda = BoxCox.lambda(y) # -0.0086589421
y_transformed = BoxCox(y, lambda)
plot(y_transformed)
lines(y_transformed) 
# Variance now stationary over time but trend component exist

# Seasonal Differencing of Period 12
y_D12 = diff(y_transformed, lag=12) # D = 1, s = 12
plot(y_D12)
lines(y_D12)
# Stationary mean, see acf and pacf for p, q, P, and Q values


# Unit root test further implies that y_D12 is stationary!
adf.test(y_D12) # pvalue=0.02 reject h0, stationary

# kpss.test(y)
# kpss.test(y_transformed)
kpss.test(y_D12) # pvalue =0.1 fair to reject h0, mean stationary



acf(y_D12, lag.max=50) # q = 5, Q = 2
pacf(y_D12, lag.max=50) # p = 10, P = 2

# Four combinations of SARIMA models
fit1 = arima(y_transformed, order=c(0,0,5), seasonal=list(order=c(0,1,2), period=12))
tsdiag(fit1) # Not adequate based on p values


# -------------FIT1 IMPROVEMENT BEGINS HERE ----------
# check residuals for improvement
acf(fit1$residuals, lag=60) # suggest residuals is q=1, Q =2 
pacf(fit1$residuals, lag=60) # suggest residuals is p=4, P=2

# try the new 4 improvements
# Fit1 Improvement 1
fit1_1 = arima(y_transformed, order=c(0,0,6), seasonal=list(order=c(0,1,4), period=12))
tsdiag(fit1_1) # still fail
acf(fit1_1$residuals, lag = 100) # residuals q is ma(2) Q is 0
pacf(fit1_1$residuals, lag = 100) # residuals p is ar(9) P is 0

# Fit1 improvement 1's improvement 1
fit1_1_1 = arima(y_transformed, order=c(0,0,8), seasonal=list(order=c(0,1,4), period=12))
tsdiag(fit1_1_1) #fail again
acf(fit1_1_1$residuals, lag = 100) # residuals q is ma(9) Q is 0
pacf(fit1_1_1$residuals, lag = 100) # residuals p is ar(9) P is 4

# Fit1 improvement 1's improvement 1 improvement 1
fit1_1_1_1 = arima(y_transformed, order=c(0,0,17), seasonal=list(order=c(4,1,4), period=12))
tsdiag(fit1_1_1_1) #pass finally
fit1_1_1_1 #AIC -499

# Fit1 improvement 1's improvement 1 improvement 2
fit1_1_1_2 = arima(y_transformed, order=c(9,0,8), seasonal=list(order=c(4,1,4), period=12))
tsdiag(fit1_1_1_2) #pass finally
fit1_1_1_2 #AIC -511.74

# Fit1 improvement 1's improvement 1 improvement 3
fit1_1_1_3 = arima(y_transformed, order=c(0,0,17), seasonal=list(order=c(0,1,4), period=12))
tsdiag(fit1_1_1_3) #pass finally
fit1_1_1_3 #AIC -490.55

# Fit1 improvement 1's improvement 1 improvement 4
fit1_1_1_4 = arima(y_transformed, order=c(9,0,8), seasonal=list(order=c(0,1,4), period=12))
tsdiag(fit1_1_1_4) #pass finally
fit1_1_1_4 #AIC -518.11 



# -

# Fit1 improvement 1's improvement 2
fit1_1_2 = arima(y_transformed, order=c(9,0,6), seasonal=list(order=c(0,1,4), period=12))
tsdiag(fit1_1_2) #pass finally
fit1_1_2 # -521.38


# Fit1 Improvement 2
fit1_2 = arima(y_transformed, order=c(4,0,5), seasonal=list(order=c(2,1,2), period=12))
tsdiag(fit1_2) # improvement passed adequate test !
fit1_2 # -519.37 very high!


# Fit1 Improvement 3
fit1_3 = arima(y_transformed, order=c(0,0,6), seasonal=list(order=c(2,1,2), period=12))
tsdiag(fit1_3) # imporvement attempt failed
acf(fit1_3$residuals, lag=30) #residuals q is ma(9), Q=0
pacf(fit1_3$residuals, lag=30) #residuals p is ar(9), P=0

# Fit 1 improvement 3's improvement 1
fit1_3_1 = arima(y_transformed, order=c(0,0,15), seasonal=list(order=c(2,1,2), period=12))
tsdiag(fit1_3_1) # pass finally
fit1_3_1 # AIC -484.51

# Fit 1 improvement 3's improvement 2
fit1_3_2 = arima(y_transformed, order=c(9,0,6), seasonal=list(order=c(2,1,2), period=12))
tsdiag(fit1_3_2) # pass finally also
fit1_3_2 # AIC -520.96 Best!


# Fit1 Improvement 4
fit1_4 = arima(y_transformed, order=c(4,0,5), seasonal=list(order=c(0,1,4), period=12))
tsdiag(fit1_4) # Pass also
fit1_4 # AIC -520

# ----------- FIT1 IMPROVEMENT ENDS HERE -------------------




fit2 = arima(y_transformed, order=c(0,0,5), seasonal=list(order=c(2,1,0), period=12))
tsdiag(fit2) # Not adequate based on p values


# ----------------- FIT2 IMPORVEMENT BEGINS HERE --------------------
# check residuals for improvement, we have 4 competiting models
acf(fit2$residuals, lag=60) # suggest residuals is q=1, Q=2
pacf(fit2$residuals, lag=60) # suggest residuals is p=4, P=2


# 2 new competiting models after checking arima(0,0,6,2,1,2)
fit2_1 = arima(y_transformed, order=c(0,0,6), seasonal=list(order=c(2,1,2), period=12))
tsdiag(fit2_1) # FAILED again. CONTINUE IMPROVEMENT FROM HERE
acf(fit2_1$residuals, lag=48) # p=9, P=0
pacf(fit2_1$residuals, lag=48) # q=9, Q=0 

# arima (9 0 6 2 1 2) already checked, as fit1_3_2
# but in this derivation we end up being implied to check arima(9 0 6 2 1 2) again
fit2_1_1 = arima(y_transformed, order=c(9,0,6), seasonal=list(order=c(2,1,2), period=12))
tsdiag(fit2_1_1) # good
fit2_1_1 # -520.96

fit2_1_2 = arima(y_transformed, order=c(0,0,15), seasonal=list(order=c(2,1,2), period=12))
tsdiag(fit2_1_2) # good
fit2_1_2 # -484.51


fit2_2 = arima(y_transformed, order=c(4,0,5), seasonal=list(order=c(4,1,0), period=12))
tsdiag(fit2_2) # good
fit2_2 # -518.46


fit2_3 = arima(y_transformed, order=c(0,0,6), seasonal=list(order=c(4,1,0), period=12))
tsdiag(fit2_3) # good
fit2_3 # -425.38


fit2_4 = arima(y_transformed, order=c(4,0,5), seasonal=list(order=c(2,1,2), period=12))
tsdiag(fit2_4) # good
fit2_4 # -519.37

# ----------------- FIT2 IMPORVEMENT ENDS HERE --------------------

fit3 = arima(y_transformed, order=c(10,0,0), seasonal=list(order=c(2,1,0), period=12))
tsdiag(fit3) # Adequate based on p values and residuals acf cut off after lag 0

fit4 = arima(y_transformed, order=c(10,0,0), seasonal=list(order=c(0,1,2), period=12))
tsdiag(fit4) # Adequate based on p values and residuals acf cut off after lag 0

# Show statistics of fit3 and fit4
fit3 # aic = -513.12
fit4 # aic = -517.68
# fit4 produced better aic, thus selected

# Plot y_transformed
plot(1:204, y_transformed, xlim=c(0,240), ylim=c(0,4), main='Drug Sales (Transformed)', xlab='Time', ylab='Drug Sales (Transformed)')
lines(y_transformed, type='l')

# Predict 12 months aheaed and plot on y_transformed plot
forecast_y_transformed = predict(fit1_1_2, n.ahead=36)
lines(205:240, forecast_y_transformed$pred, type="o", col="red")

# inverse transform 
forecast_y = InvBoxCox(forecast_y_transformed$pred, lambda)
forecast_ub = InvBoxCox(forecast_y_transformed$pred + 1.96*forecast_y_transformed$se, lambda) # upper bound
forecast_lb = InvBoxCox(forecast_y_transformed$pred - 1.96*forecast_y_transformed$se, lambda) # lower bound

# Plot original y data and forecast + upperbound and lower bound of forecast predictions
plot(1:204, y, xlim=c(0,240), ylim=c(0,55), main="Drug Sales", xlab='Month', ylab='Drug Sales')
lines(y, type='l')
lines(205:240, forecast_y, type='l', col='red')
lines(205:240, forecast_ub, type='l', col='blue')
lines(205:240, forecast_lb, type='l', col='blue')
legend("topleft", c("Drug Data", "3-years Ahead Forecast", "Forecast Upper & Lower Bounds"), col=c('black','red','blue'), lty=c(1,1,1))

#Holt-Winters' Trend and Seasonality Model for prediction

install.packages("dplyr")
install.packages("zoo")
install.packages("ggplot2")
library(dplyr)
library(zoo)
library(ggplot2)

y2 = read.delim("drug.txt", sep = ",", header = TRUE)
global.freq <<- 12 # seasonality
y2 = y2 %>% mutate(date = ymd(date))
global.start <<- ymd(as.Date(y2$date[[1]]))
global.end <<- ymd(as.Date(y2$date[[nrow(y2)]]))
### Convert to time series data
y2_all = ts(y2$value, start=c(year(global.start), month(global.start)),
            end=c(year(global.end), month(global.end)), frequency=global.freq)

plot(y2_all, main='Time Series: Drug') 
no_of_train = round(0.8 * length(y2_all))
y2_train = head(y2_all, no_of_train)
y2_test = tail(y2_all, round(length(y2_all) - no_of_train))
y2 = y2$value
lambda3 = BoxCox.lambda(y2) # -0.0086589421

lambda2 = BoxCox.lambda(y2_all) # 0.1313326
# lambda2 = 0
# lambda3 = 0

# lambda2 = BoxCox.lambda(y2_all) # 0.1313326
# lambda2 = 0
y2_all
y2

lambda2 = BoxCox.lambda(y2_all) # 0.1313326
lambda3 = 0


y2_all_t = BoxCox(y2_all,lambda3) 
y2_train_t = BoxCox(y2_train,lambda3) 
y2_test_t = BoxCox(y2_test,lambda3) 

fit.hw.add = hw(y2_train_t, seasonal = "additive")
fit.hw.mul = hw(y2_train_t, seasonal = "multiplicative")

pred.hw.add = InvBoxCox(forecast(fit.hw.add$mean, h=length(y2_test_t)+12)$mean, lambda = lambda3) # +12 so that the final date of prediction is 3 yrs ahead of data provided
pred.hw.mul = InvBoxCox(forecast(fit.hw.mul$mean, h=length(y2_test_t)+12)$mean, lambda = lambda3) # +12 so that the final date of prediction is 3 yrs ahead of data provided

global.pred.start = start(pred.hw.add) %>% print()
global.pred.end = end(pred.hw.add) %>% print()

acc.hw.add = accuracy(pred.hw.add, y2_all_t) %>% print()
acc.hw.mul = accuracy(pred.hw.mul, y2_all_t) %>% print()

#Dataframe for all data
df1 = data.frame(value=as.matrix(y2_all), date=as.Date(as.yearmon(time(y2_all))), model = "Base")
args = list(pred.hw.add[18:53], "HW's ADD", pred.hw.mul[18:53], "HW's MUL")
if (!length(args) == 0) {
  for(i in seq(1, length(args), by=2)) {
    name = toString(args[[i+1]])
    # forecast = ts(args[[i]], start=c(global.pred.start), end=c(global.pred.end), frequency=global.freq)
    forecast = ts(args[[i]], start=c(2008, 7), end=c(2011, 6), frequency=global.freq)
    df <- data.frame(value=as.matrix(forecast), date=as.Date(as.yearmon(time(forecast))), model = name)
    df1 <- rbind(df1, df)
  }
}

# Original sky plot
ggplot(df1, aes(date, value, colour = model)) +geom_line()+geom_vline(aes(xintercept = as.numeric(date[length(y2_all)])),
                                                                      linetype = "longdash", color = "black")+
  ylab("")+
  xlab("")+
  ggtitle("Time Series: Holt-Winter Model forecast")


#only show the forecast in last 3 year with last 50 values from train dataset
df2 = data.frame()
df2 = data.frame(value=as.matrix(tail(y2_all,50)), date=as.Date(as.yearmon(time(tail(y2_all,50)))), model = "Base")
args = list(pred.hw.add, "HW's ADD", pred.hw.mul, "HW's MUL")
if (!length(args) == 0) {
  for(i in seq(1, length(args), by=2)) {
    name = toString(args[[i+1]])
    forecast = ts(args[[i]], start=c(global.pred.start), end=c(global.pred.end), frequency=global.freq)
    df <- data.frame(value=as.matrix(forecast), date=as.Date(as.yearmon(time(forecast))), model = name)
    df2 <- rbind(df2,df)
  }
}


ggplot(df2, aes(date, value, colour = model)) +geom_line()+geom_vline(aes(xintercept = as.numeric(date[length(y2_test_t)])),
                                                                      linetype = "longdash", color = "black")+
  ylab("")+
  xlab("")+
  ggtitle("Time Series: Holt-Winter Model forecast-3year")



# plot comparison of forecast
plot(1:36, forecast_y, xlim=c(0,36), ylim=c(15,50), main='Forecast Comparison', xlab='Months Ahead', ylab='Drug Sales', type='l', col='red')
lines(1:36, pred.hw.add[18:53], type='l', col='green')
lines(1:36, pred.hw.mul[18:53], type='l', col='blue')
legend("topleft", c("SARIMA((9,0,6),(0,1,4),12)", "HW Add", "HW Mul"), col=c('red','green', 'blue'), lty=c(1,1,1))


plot(1:204, y, xlim=c(0,240), ylim=c(0,50), main='Forecast Comparison', xlab='Months Ahead', ylab='Drug Sales', type='l', col='black')
lines(205:240, forecast_y, type='l', col='red')
lines(205:240, pred.hw.add[18:53], type='l', col='green')
lines(205:240, pred.hw.mul[18:53], type='l', col='blue')
legend("topleft", c("Drug Sales Data", "SARIMA((9,0,6),(0,1,4),12)", "HW Add", "HW Mul"), col=c('black', 'red','green', 'blue'), lty=c(1,1,1))
