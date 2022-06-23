#' Detects when the laser blasts through the target sample in a laser ablation mass spectrometry time resolved analysis.
#'
#' This function imports a data frame containing a single time resolved laser ablation mass spectrometry analysis of a foraminifera (or other carbonate shell).
#' It assumes that the first row contains the signal of the target sample and that background correction has already been applied.
#' \cr
#' Column names referencing the time stamps and target signal are specified as function arguments.
#' Then the signal column is evaluated by smoothing the values using a moving average, then scaled between 0-1 and the rate of change over a number of observations.
#' \cr
#' The function identifies the maximum rate of signal change, locates the corresponding time stamp and then subtracts the time it took for the laser to blast through the target.
#' The result (endTime) is the time stamp of the last observation while the laser is still focussed the desired target.
#'
#' @param df A data frame containing a single time resolved analysis, with a column referencing time and another with the corresponding measured counts data.
#' @param dt An integer that controls the number of observations (rows) are used in calculating a rolling lagged difference in 44Ca signal. Using a lower value for a faster blast through of chamber wall can improve end point sensitivity. Default = 10.
#' @param smoothing Controls the length of the moving average filter to apply over the dt period. Default = 5.
#' @param timeCol The column title in the data frame identifying the time stamp of the time resolved analysis. Default = "Time".
#' @param signalCol The column title in the data frame identifying the numerical data used to identify the laser ablation blast through endpoint. This is often 44Ca but could be any column of numerical data that you want to detect the endpoint. Default = "Ca44".
#' @param profile Logical. A visualisation of the endpoint detection mechanism as a ggplot2 object. To make this profile plot, set the argument to TRUE, otherwise set it to FALSE. Specifying to not make a plot can save a substantial amount of time. Default = "FALSE".
#' @param timeUnits The units that the time resolved analysis is measured in. This is the units of the timeCol. This argument is a string and is only necessary if the argument profile = "TRUE". Default = "seconds".
#'
#' @return The function returns a data frame containing the columns:
#' * dfReturn$df contains a data frame with only the observations between the first data frame row and the endTime.
#' * dfReturn$startTime contains the earliest time in your TRA as a numerical value.
#' * df$Return$endTime contains the last time step while the laser is still focussed the desired target in your TRA before the as a numerical value.
#' * df$Return$profile contains a visualisation of your TRA identifying where the laser ablated through the carbonate shell as a ggplot object. This is only available if a profile was generated using profile = "TRUE".
#'
#'
#' @examples
#' endPoint(df = foram72shot3, dt = 10, smoothing = 5, timeCol = "Time",
#'  signalCol = "Ca44", profile = "TRUE",  timeUnits = "seconds")
#' endPoint(df = foram166shot7, dt = 8, smoothing = 7, timeCol = "Time",
#'  signalCol = "Ca44", profile = "FALSE",  timeUnits = "seconds")
#' endPoint(df = foram174shot4, dt = 10, smoothing = 5, timeCol = "Time",
#'  signalCol = "Ca43", profile = "TRUE",  timeUnits = "seconds")
#' endPoint(df = coral6, dt = 10, smoothing = 5, timeCol = "Time",
#'  signalCol = "Sr86", profile = "FALSE",  timeUnits = "milliseconds")

#' @export
#'
#' @importFrom rlang .data

endPoint <- function(df, dt = 10, smoothing = 5, timeCol = "Time", signalCol = "Ca44", profile = "FLASE",  timeUnits = "seconds"){

  # scanRate is the frequency of signal detection. The number of rows of data per second in your data frame.
  scanRate <- 1 / (df$Time[2] - df$Time[1])

  # smooth the signal before detecting the rate of change.
  # A large value of order causes more smoothing
  df$Ca44sma <- smooth::sma(df[, grep(signalCol, names(df))[1]], order = smoothing)$fitted %>% as.numeric()

  # This gives change as absolute
  Ca44dydt <- diff(df$Ca44sma, lag = dt) / diff(df[, grep(timeCol, names(df))[1]], lag = dt)
  Ca44dydt <- c(rep(0, times = dt), Ca44dydt)
  Ca44dydt <- (Ca44dydt - min(Ca44dydt)) / (max(Ca44dydt) - min(Ca44dydt)) #scaled to 0-1
  df$Ca44dydt <- Ca44dydt

  #The scaled values (between 0 - 1) of 44Ca are calculated here for plotting on a graph
  #for rate of change against time ($Time) comparison
  df$Ca44Scaled <- (df$Ca44sma - min(df$Ca44sma))/(max(df$Ca44sma) - min(df$Ca44sma))

  # the number of seconds to remove after end point detection (endTime)
  tailSeconds <- (dt / scanRate) # the number of seconds to remove after end point detection (endTime)
  tailSeconds <- round(tailSeconds, digits = 3) #not working with pipes

  # Identify the time ($Time) in seconds when the maximum rate of change in signalCol occurs.
  # The subtraction of tailSeconds is in place to backtrack the zero dydt values at the start. This results in the endTime being before the sharp signal cutoff.
  # Notice that dt is in observations, but endTime in in seconds.
  endTime <- subset(df[, grep(timeCol, names(df))[1]], df$Ca44dydt == min(df$Ca44dydt)) - tailSeconds

  # A consequence of legacy code
  startTime <- min(df[, grep(timeCol, names(df))[1]])
  maxTime <- max(df[, grep(timeCol, names(df))[1]])

  # Make a data frame to append the outputs from this function to.
  # This data frame will be returned outside this function
  dfReturn <- NULL

  dfReturn$df <- df # Your data frame with only the rows that are between your startTime and endTime
  dfReturn$exceededThresholdSubset <- NULL #not currently used
  dfReturn$startTime <- startTime # The first time step in your analysis after this function
  dfReturn$endTime <- endTime # The last time step in your analysis after this function has been applied.

  if(profile == "TRUE"){
    dfReturn$profile <- ggplot2::ggplot(df, ggplot2::aes(x=df[, grep(timeCol, names(df))[1]])) +
      ggplot2::annotate("rect", xmin = startTime - 10, xmax = startTime, ymin = -Inf, ymax = Inf, fill ="red", alpha = 0.5)+ #change -10 to a percentage
      ggplot2::annotate("rect", xmin = endTime, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.5)+
      ggplot2::geom_point(ggplot2::aes(y=df$Ca44Scaled, colour = "Signal")) +
      ggplot2::geom_line(ggplot2::aes(y=df$Ca44Scaled, colour = "Signal")) +
      ggplot2::geom_line(ggplot2::aes(y=df$Ca44dydt, colour = "dydt")) +
      ggplot2::geom_vline(xintercept = endTime, colour = "purple")+
      ggplot2::geom_label(x = startTime - 10 , y = stats::median(abs(scale(df[, grep(signalCol, names(df))[1]], center = TRUE))), label = paste("TRA started at", startTime, timeUnits), size = 2, hjust = "left")+
      ggplot2::geom_label(x = endTime, y = mean(abs(scale(df[, grep(signalCol, names(df))[1]], center = TRUE))), label = paste("endTime \n (largest signal change - (dt/scanRate)) \n at", endTime, timeUnits), size = 2)+
      ggplot2::labs(y = paste("Scaled", signalCol, "signal and rate of change"),
           x = paste("Time elapsed in", timeUnits),
           subtitle = paste("With a smoothing of", smoothing, "observations, \n a dt of", dt, "observations \n and dt/scanRate of", tailSeconds, timeUnits)
      )+
      ggplot2::scale_x_continuous(limits = c(startTime-10, maxTime),
                         breaks = scales::pretty_breaks(n = 10))+
      ggplot2::scale_colour_manual("",
                          breaks = c("Signal", "dydt"),
                          values = c("black", "blue")) +
      ggplot2::theme_bw()

    # Data points located within the red shading are removed in the returned data frame as these are beyond

    # The user can add a title to this plot after running this function by adding the following lines of code
    # dfReturn$profile + labs(title = "Title text")

  }

  # Return data frame to global environment
  #.GlobalEnv$dfReturn <- dfReturn

  return(dfReturn)

}
