FROM rocker/shiny

LABEL image.author.name="MartonHorvath98"
LABEL image.author.email="marton.hrvth.98@gmail.com"

# Install MySQL client libs
RUN apt-get update && apt-get install -y default-libmysqlclient-dev

# Make a directory in the container
RUN mkdir /home/shiny-app
RUN mkdir -p /home/shiny-app/meta
RUN mkdir -p /home/shiny-app/www

ENV DB_HOST=host.docker.internal
ENV DB_USER=root
ENV DB_NAME=InteractomeDB
# ENV DB_PASSWORD=!!!pass-when-running

# Install R dependencies
RUN R -e "install.packages(c('shiny', 'DBI', 'RMySQL', 'DT', 'dplyr', 'openxlsx'))"

# Copy the Shiny app code
COPY app.R /home/shiny-app/app.R

# Copy meta folder and www
COPY meta/SURFME_v2023.xlsx /home/shiny-app/meta/SURFME_v2023.xlsx
COPY www/nsfw-logo.png /home/shiny-app/www/nsfw-logo.png

WORKDIR /home/shiny-app

# Expose the application port
EXPOSE 8180

# Run the R Shiny app
CMD ["R", "-e", "shiny::runApp('/home/shiny-app/app.R', port=8180, host='0.0.0.0')"]