
##############################################################################
#                                                                            #
#    U.S. Monthly Wind Power Production Forecasting                          #
#                                                                            #
############################################################################## 

################################## Problem (a) #################################
library(zoo)
library(ggplot2)
library(ggpubr)
library(forecast)

data <- read.csv("wind-power-production-us.csv") 
dim(data)

wind_power_production <- ts(data$wind_united_states, start = c(2001, 1), freq = 12)
wind_power_production_log <- ts(log(data$wind_united_states), start = c(2001, 1), freq = 12)

fig_1 <- autoplot.zoo(wind_power_production) + 
  xlab("Year")+
  ylab("Wind Power production")
ggtitle("Wind Power production in USA in 2001 to 2023")
xlim(2001, 2023+2/12)
plot(fig_1)

fig_2 <- autoplot.zoo(wind_power_production_log) + 
  xlab("Year")+
  ylab("Log-transformed Wind Power production")
ggtitle("Log-transformed Wind Power production in USA in 2001 to 2023")
xlim(2001, 2023+2/12)
plot(fig_2)

plot_wind <- ggarrange(fig_1, fig_2, nrow = 2, ncol = 1, align = "v")
plot_wind




################################## Problem (b) #################################

library(smoots)
library(forecast)


ts_zoo <- zoo(wind_power_production_log)
rate <- as.ts(ts_zoo)

est <- msmooth(rate)
bwidth <- est$b0
bwidth


plot_rate <- autoplot.zoo(rate) +
  xlab("Year") +
  ylab("Log-transformed data of Wind Power production rate") +
  ggtitle("Log transformed data of Wind Power production")
plot(plot_rate)  

trend <- fitted(est)
df <- data.frame(
  t = time(rate),
  trend = trend
)

plot_trend <- plot_rate +
  geom_line(data = df, aes(x = t, y = trend), color = "red", linewidth = 0.8) +
  ggtitle("The observed series (black) together with the estimated local line trend")

plot_trend

################################## Problem (c) #################################

res <- resid(est)
plot_residuals <- autoplot.zoo(res) +
  xlab("Year") +
  ylab("Residuals Values")+
ggtitle("The residuals series")+
  geom_hline(yintercept = 0, color ="red")
plot_residuals

acf <- ggAcf(as.numeric(res))+
  ggtitle("Correlogram of the deterended series")

acf



################################## Problem (d) #################################


p_max <- q_max <- 2
p <- 0:p_max
q <- 0:q_max
bic <- matrix(NA, nrow = p_max + 1, ncol = q_max + 1)
rownames(bic) <- paste0("p=", p)
colnames(bic) <-  paste0("q=", q)
n <- length(res)

for (p0 in p) {
  for (q0 in q) {
    arma <- arima(res, order = c(p0, 0, q0), include.mean = FALSE)
    bic[(p0 + 1), (q0 + 1)] <- AIC(arma, k = log(n))
  }
  
}
bic

pq_opt <- unname(which(bic == min(bic), arr.ind = TRUE) - 1)
p_opt <- pq_opt[[1]]
q_opt <- pq_opt[[2]]
p_opt
q_opt


arma_opt <- arima(res,order = c(p_opt, 0, q_opt), include.mean = FALSE)
arma_opt


################################## Problem (e) #################################


n <- length(wind_power_production)    # Number of observations
# Number of test observations
n_te <- trunc(0.1 * n)
# Number of training observations
n_tr <- n - n_te 

# Absolute bandwidth
b_abs <- trunc(est$b0 * n_tr + 0.5)
# Readjusted relative bandwidth for complete data
b_n <- b_abs / n




est_n <- gsmooth(wind_power_production_log, b = b_n)
arma_n <- Arima(est_n$res, order = c(p_opt, 0, q_opt), include.mean = FALSE)
arma_n

arma_res <- arma_n$residuals




library(tseries)

jarque.bera.test(res)

# A bootstrap requires setting a seed for reproducibility

set.seed(123)

#Forecast for the parametic part
fc_para <- forecast(arma_n, h = 30, bootstrap = TRUE, level = 95)


# Forecast for the non-parametic part
fc_trend <- tail(est_n$ye, 1) + (1:30) * diff(tail(est_n$ye, 2))

# Forecast according to the complete model
fc_point <- exp(fc_para$mean + fc_trend)
fc_low <- exp(fc_para$lower + fc_trend)
fc_up <- exp(fc_para$upper + fc_trend)

# data with the results

df <- data.frame(
   Year = c(time(fc_point)),
            fc_point = c(fc_point),
            fc_low = c(fc_low),
            fc_up = c(fc_up)
          )
df

# Plot the original series with point and interval forecast at the end

autoplot.zoo(wind_power_production) +
  geom_ribbon(
    data = df,
    aes(x= Year, ymin = fc_low, ymax = fc_up),
    fill = alpha("red", 0.5)
  ) +
  geom_line(data = df, aes(x = Year, y = fc_point), color = "blue") +
  ggtitle("Forecast of the Wind Power production in USA") +
  xlab("Year") +
  ylab("Wind Power production")


# Forecast for the whole Semi-Arima model

fc_mc <- exp(modelCast(est_n, p = p_opt, q = q_opt, h = 30,
                       method = "boot"))
fc_mc

################################ END ################################




