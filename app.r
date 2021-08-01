library('shiny')
library('tidyverse')
library('shinythemes')
library('shinyWidgets')
library('plotly')
library('maps')
library('coda')

setwd('D:\\Advanced Statistics for Physics Analysis lecture material\\esame\\webapp\\source\\mio2')

if (!file.exists('./mcmc_data/mcmc_data.RData')) {
    source('mcmc_data_generator.r')
}
load('./mcmc_data/mcmc_data.RData')
load('./mcmc_data/mcmc_data_by_age.RData')
# eff_list <- ls() # this lazy option is dangerous if the terminal session isn't restarted every time...
eff_list <- list(eff.pfizer.old = eff.pfizer.old, eff.pfizer.new = eff.pfizer.new, eff.moderna = eff.moderna, eff.astrazeneca = eff.astrazeneca, eff.janssen = eff.janssen)
param_list <- list(
                   pv.pfizer.old = pv.pfizer.old, pp.pfizer.old = pp.pfizer.old, 
                   pv.pfizer.new = pv.pfizer.new, pp.pfizer.new = pp.pfizer.new, 
                   pv.moderna = pv.moderna, pp.moderna = pp.moderna, 
                   pv.astrazeneca = pv.astrazeneca, pp.astrazeneca = pp.astrazeneca, 
                   pv.janssen = pv.janssen, pp.janssen = pp.janssen
)

source('mcmc_utilities.R')

UPDATE_DATA <- TRUE # set to TRUE to redownload the owid csv
if (UPDATE_DATA) {
    write.csv(read.csv('https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/vaccinations/vaccinations.csv'), 'vaccinations.csv')
}
df <- read.csv('vaccinations.csv')

df$date <- as.Date(df$date, '%Y-%m-%d') # non fondamentale con ggplot ma silenzia un warning e va più veloce

MIN_DATE <- min(df$date)
MAX_DATE <- max(df$date)

print(paste('most recent date in the dataset:', MAX_DATE))
print(paste('today is', lubridate::today()))

plot_vaccinations_stats <- function(data, loc, col, dates) {
    tmp <- data |> drop_na() |> filter(location %in% loc) |> filter(date > dates[1] & date < dates[2])# |> select(col)
    g <- ggplot(tmp, aes_string(x = 'date', y = col, color = 'location')) + geom_line() + geom_point() + xlab('date') + ylab(str_replace(col, '_', ' '))
    g
}

plot_world_map <- function(world_data) {
    a = world_data |> select(location, date, people_vaccinated_per_hundred) |> drop_na() |> group_by(location) |>
                summarise(maxDate = max(date), vacc = people_vaccinated_per_hundred[which.max(date)]) |> rename(region = location)
    world = map_data('world') |> mutate(region = replace(region, region == 'USA', 'United States')) |> mutate(region = replace(region, region == 'UK', 'United Kingdom'))  
    # it is useful in order to plot with ggplot. It is just a map of the world with column about longitude and latitude, regions and subregions
    world_map = left_join(world, a, by = 'region') 
    # left_join uses a kinda of groupby operation with the "by" option that allows us to pair informations from world and a by matching regions
    #options(repr.plot.width=20, repr.plot.height=10)
    # changing plot options 
    g <- ggplot(data = world_map) + 
         geom_polygon(aes(x = long, y = lat, fill = vacc, group = group), color = "white") +
         coord_fixed(1.3) + ggtitle(paste('Percentage of vaccinated people per nation (at least one dose), last updated:', MAX_DATE)) +
         theme(plot.title = element_text(color="black", size=16, face="bold", hjust = 0.5),
         legend.text = element_text(size=16)) + labs(fill = 'Percentage')
    g
}

pfizer.df <- read.csv('./vaccine_trials_data/pfizer_data.csv')
moderna.df <- read.csv('./vaccine_trials_data/moderna_data.csv')
astrazeneca.df <- read.csv('./vaccine_trials_data/astrazeneca_data.csv')
janssen.df <- read.csv('./vaccine_trials_data/janssen_data.csv')

pfizer.age.df <- read.csv('./age_data/pfizer_age.csv')
moderna.age.df <- read.csv('./age_data/moderna_age.csv')
janssen.age.df <- read.csv('./age_data/janssen_age.csv')
age_df_list <- list('pfizer' = pfizer.age.df, 'moderna' = moderna.age.df, 'janssen' = janssen.age.df)
age_fun_list <- list('pfizer' = pfizer.age.fun, 'moderna' = moderna.age.fun, 'janssen' = janssen.age.fun)

# conjugate prior comparison variables
eff_list_named <- list('pfizer_7d' = eff.pfizer.old, 'pfizer_6m' = eff.pfizer.new, 'moderna' = eff.moderna, 'astrazeneca' = eff.astrazeneca, 'janssen' = eff.janssen)
source('create_dataframes_for_cp_comp.r')
pv_data_named <- list('pfizer_7d' = pfizer.data.pv.old, 'pfizer_6m' = pfizer.data.pv.new, 'moderna' = moderna.data.pv, 'astrazeneca' = astrazeneca.data.pv, 'janssen' = janssen.data.pv)
pp_data_named <- list('pfizer_7d' = pfizer.data.pp.old, 'pfizer_6m' = pfizer.data.pp.new, 'moderna' = moderna.data.pp, 'astrazeneca' = astrazeneca.data.pp, 'janssen' = janssen.data.pp)


ui <- fluidPage(
    navbarPage(
        'covid-19 vaccinations tracker',
        theme = shinytheme('paper'), # paper o flatly, darkly o journal
        tabPanel(
            'map',
            plotlyOutput('map_plot')
        ),
        tabPanel(
            'vaccination stats plots',
            sidebarLayout(
                sidebarPanel(
                    pickerInput(
                        'location_select',
                        'Country/Region:',
                        choices = unique(df$location),
                        selected = c('Italy', 'United States', 'France', 'Germany', 'United Kingdom'),
                        multiple = TRUE
                    ),
                    pickerInput(
                        'col_select',
                        'Variable:',
                        choices = (names(df))[-(1:4)],
                        multiple = FALSE
                    ),
                    dateRangeInput(
                        'plot_dates',
                        'Select min./max. plot dates:',
                        start = MIN_DATE,
                        end = MAX_DATE,
                        min = MIN_DATE,
                        max = MAX_DATE,
                        language = 'en'# ,
                        # format = "yyyy-mm-dd" # default
                    )
                ),
                mainPanel(
                    plotlyOutput('plotly_plot')
                    # tabsetPanel(
                    #     tabPanel('Cumulative', plotlyOutput('plots_cumulative')),
                    #     tabPanel('New', plotlyOutput('plots_new')),
                    #     tabPanel('Cumulative (log10)', plotlyOutput('plots_cumulative_log'))
                    # )
                )
            )
        ),
        tabPanel(
            'covid vaccinations data',
            checkboxInput('include_na_in_data', 'Include NA values', value = FALSE),
            downloadButton('downloadData', 'Download data as CSV'),
            dataTableOutput('covid_data'),
            #'source: https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/vaccinations/vaccinations.csv'
        ),
        tabPanel(
            'raw vaccinations data',
            tableOutput('raw_covid_data')
        ),
        tabPanel(
            'covid vaccine trials results',
            tabsetPanel(
                tabPanel(
                    'Pfizer',
                    tableOutput('pfizer_data'),
                    '"tot_7d" and "tot_6m" refer to data collected 7 days/6 months after the second dose irrespective of patient age',
                    'Sources:',
                    tags$a(href = 'https://www.fda.gov/media/144413/download', 'FDA factsheet'),
                    tags$a(href = 'https://clinicaldata.ema.europa.eu/web/cdp/home', 'EMA clincal overview'),
                    tags$a(href = 'https://www.pfizer.com/news/press-release/press-release-detail/pfizer-and-biontech-confirm-high-efficacy-and-no-serious', 'Pfizer website (1)'),
                    tags$a(href = 'https://www.pfizer.com/news/press-release/press-release-detail/pfizer-biontech-announce-positive-topline-results-pivotal', 'Pfizer website (2)')
                ),
                tabPanel(
                    'Moderna',
                    tableOutput('moderna_data'),
                    'Sources:',
                    tags$a(href = 'https://www.fda.gov/media/144637/download', 'FDA factsheet'),
                    tags$a(href = 'https://investors.modernatx.com/news-releases/news-release-details/moderna-announces-teencove-study-its-covid-19-vaccine', 'Moderna website')
                ),
                tabPanel(
                    'Astrazeneca',
                    tableOutput('astrazeneca_data'),
                    'Sources:',
                    tags$a(href = 'https://www.thelancet.com/action/showPdf?pii=S0140-6736%2821%2900432-3', 'Oxford & collaborators study')
                ),
                tabPanel(
                    'Janssen',
                    tableOutput('janssen_data'),
                    'Sources:',
                    tags$a(href = 'https://www.fda.gov/media/146304/download', 'FDA factsheet'),
                )
            )
        ),
        tabPanel(
            'vaccines parameters (MCMC)',
            selectInput('param_name_plot', label = 'Select vaccine parameters to infer via JAGS', choices = names(param_list)),
            plotOutput('param_plot')
        ),
        tabPanel(
            'VE (MCMC)',
            selectInput('eff_name_plot', label = 'Select vaccine effectiveness to infer via JAGS', choices = names(eff_list)),
            plotOutput('eff_plot')
        ),
        tabPanel(
            'VE comparison',
            plotOutput('eff_comparator_plot')
        ),
        tabPanel(
            'VE by age',
            selectInput('age_select', label = 'Select vaccine:', choices = c('pfizer', 'moderna', 'janssen')), # age_df_list
            tableOutput('age_table'),
            plotOutput('age_plot')
        ),
        tabPanel(
            'MCMC vs CP',
            selectInput('vaccine_cp_select', label = 'Select vaccine:', choices = c('pfizer_7d', 'pfizer_6m', 'moderna', 'janssen')),
            plotOutput('vaccine_cp_plot')
        ),
        tabPanel(
            'beta vs uniform priors',
            selectInput('vaccine_priors_select', label = 'Select vaccine:', choices = c('pfizer_7d', 'pfizer_6m', 'moderna', 'janssen')),
            plotOutput('vaccine_priors_plot')
        )    
    )
)


server <- function(input, output, session) {
    output$covid_data <- renderDataTable({
        if (input$include_na_in_data) {
            df
        } else {
            df |> drop_na()
        } 
    })
    output$raw_covid_data <- renderTable(df)
    output$downloadData <- downloadHandler(
        filename = 'vaccinations.csv',
        content = \(file) write.csv(df, file)
    )
    output$plotly_plot <- renderPlotly(plot_vaccinations_stats(df, input$location_select, input$col_select, input$plot_dates))

    output$pfizer_data <- renderTable({
        pfizer.df |> rename(age = X, 'n. of positive vaccinated people' = vpos, 'n. of positive unvaccinated people (placebo)' = ppos, 'tot. n. of vac. people' = vtot, 'tot. n. of people with placebo' = ptot)
    })
    output$moderna_data <- renderTable({
        moderna.df |> rename(age = X, 'n. of positive vaccinated people' = vpos, 'n. of positive unvaccinated people (placebo)' = ppos, 'tot. n. of vac. people' = vtot, 'tot. n. of people with placebo' = ptot)
    })
    output$astrazeneca_data <- renderTable({
        astrazeneca.df |> rename(age = X, 'n. of positive vaccinated people' = vpos, 'n. of positive unvaccinated people (placebo)' = ppos, 'tot. n. of vac. people' = vtot, 'tot. n. of people with placebo' = ptot)
    })
    output$janssen_data <- renderTable({
        janssen.df |> rename(age = X, 'n. of positive vaccinated people' = vpos, 'n. of positive unvaccinated people (placebo)' = ppos, 'tot. n. of vac. people' = vtot, 'tot. n. of people with placebo' = ptot)
    })

    output$param_plot <- renderPlot(diagMCMC(codaObject = param_list[[input$param_name_plot]], parName = "theta", title = input$param_name_plot)) 
    output$eff_plot <- renderPlot(plotPost((eff_list[[input$eff_name_plot]])*100, main = input$eff_name_plot, xlab = '', hdicol = 'cornflowerblue'))

    output$map_plot <- renderPlotly(plot_world_map(df))

    output$eff_comparator_plot <- renderPlot(barcomparator(list(eff.pfizer.old, eff.pfizer.new, eff.moderna, eff.astrazeneca, eff.janssen)))

    output$age_plot <- renderPlot(age_fun_list[[input$age_select]]())
    output$age_table <- renderTable(age_df_list[[input$age_select]] |> rename(age = X))

    output$vaccine_cp_plot <- renderPlot({
        plotPost1((eff_list_named[[input$vaccine_cp_select]])*100, posterior(df.pv = pv_data_named[[input$vaccine_cp_select]], df.pp = pp_data_named[[input$vaccine_cp_select]]), 
        main = paste0("Posterior comparison (", input$vaccine_cp_select, '): MCMC vs MC + kde'), xlab = "VE", hdicol1 = 'dodgerblue', hdicol2 = "red") # , showCurve = TRUE
    })

    output$vaccine_priors_plot <- renderPlot(
        priorcomparison(
            posterior(df.pv = pv_data_named[[input$vaccine_priors_select]], df.pp = pp_data_named[[input$vaccine_priors_select]]),
            posterior.uniform(df.pv = pv_data_named[[input$vaccine_priors_select]], df.pp = pp_data_named[[input$vaccine_priors_select]]),
            xlab = 'VE',
            main = paste0('Posterior comparison (', input$vaccine_priors_select, '): beta vs uniform priors')
        )
    )

}

runApp(shinyApp(ui, server))