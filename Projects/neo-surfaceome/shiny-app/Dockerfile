FROM rocker/shiny

LABEL image.author.name="MartonHorvath98"
LABEL image.author.email="marton.hrvth.98@gmail.com"

# Make a directory in the container
RUN mkdir /home/shiny-app

# Install R dependencies
RUN R -e "install.packages(c('dplyr', 'ggplot2', 'gapminder'))"

# Copy the Shiny app code
COPY app.R /home/shiny-app/app.R

# Expose the application port
EXPOSE 8180

# Run the R Shiny app
CMD Rscript /home/shiny-app/app.R