
####### for running in local web browser ###############
import nest_asyncio
nest_asyncio.apply()
########################################################

####### for functions ##################################
import asyncio
import os
import pandas as pd
import numpy as np
#import datatable as dt
#import xlsxwriter
########################################################

################### For Shiny ##########################
import urllib.request
from pathlib import Path

import duckdb
from query import query_output_server, query_output_ui
#from query import query_output_ui
from shiny import App, reactive, ui, render
#######################################################
#for local dev:
os.chdir("C:/Users/mmari/OneDrive/Documents/GitHub/qbs_app/")

app_dir = Path(__file__).parent
db_file = app_dir / "weather.db"

def load_csv(con, csv_name, table_name):
    csv_url = f"https://raw.githubusercontent.com/rfordatascience/tidytuesday/\
        master/data/2022/2022-12-20/{csv_name}.csv"
    local_file_path = app_dir / f"{csv_name}.csv"
    urllib.request.urlretrieve(csv_url, local_file_path)
    con.sql(
        """
        f"CREATE TABLE {table_name} AS SELECT * 
        FROM read_csv_auto('{local_file_path}')"
        """
    )

if not Path.exists(db_file):
    con = duckdb.connect(str(db_file), read_only=False)
    load_csv(con, "weather_forecasts", "weather")
    load_csv(con, "cities", "cities")
    con.close()

con = duckdb.connect(str(db_file), read_only=True)

app_ui = ui.page_sidebar(
    ui.sidebar(
        ui.input_action_button("add_query", 
                               "Upload schedule file", 
                               class_="btn btn-primary"),
        ui.div("(working on final schedule file format ...)"),
        ui.input_action_button("run_matches", "Run Matching Routine"),
        ui.download_button("downloadData", 
                           "Download Matches", 
                           class_="btn btn-secondary"),
        #ui.input_action_button(
        #    "show_meta", "Dowload Results", 
        #),
        #ui.markdown(
        #    """
        #    This App is for matching QBS applicants with desired faculty 
        #    and alumni
        #    """
        #    #This app lets you explore a dataset using SQL and duckdb.
        #    #The data is stored in an on-disk [duckdb](https://duckdb.org/) 
        #    database,
        #    #which leads to extremely fast queries.
        #),
        ui.input_radio_buttons(
            "radio_group",  # Input ID
            #"Select an option:",  # Label for the radio button group
            "Show current data:",
            ["Show applicants", 
             "Show faculty", 
             "Show alumni",
             "Show times",
             "Show faculty times",
             "Show days"],  # List of choices
            ),
    ),
    ui.tags.div(
        ui.layout_column_wrap(
        #ui.layout_columns(
            #ui.tags.style(
            #    ".output_text { border: 1px solid #ccc; padding: 5px; 
            #    border-radius: 4px; background-color: #f8f8f8; 
            #    font-family: monospace; white-space: pre-wrap;}"
            #),
            #ui.output_text("radio_value"),
            ui.card(
                ui.card_header("Output selected file contents here"),
                #ui.p("This is the body."),
                ui.output_text("radio_value"),
                ui.output_text_verbatim("radio_frame"),
                #ui.p("This is still the body."),
                #ui.card_footer("This is the footer"),
                full_screen=True,
            ),
            ui.card(
                #ui.card_header("Data Frame as ", 
                #ui.tags.code("render.DataTable")),
                ui.card_header("Output final matches here"),
                ui.output_data_frame("table"),
            ),
            width=1 / 2,
            max_height="500px",
         #   query_output_ui("initial_query", remove_id="initial_query"),
          #  id="module_container",
            #col_widths = {"lg": [4, 8]},
        #),
        ),
    ),
    #title="DuckDB query explorer",
    title='''
        "QBS App - this App is for matching QBS applicants 
        with desired faculty and alumni"
        ''',
    class_="bslib-page-dashboard",
)


def server(input, output, session):
    
    
    ##################### Main inut data (.xlsx) #############################
    ##########################################################################
    ##########################################################################
    
    inputPath='input_file.xlsx'
    def get_input_file(inputPath):
        inputPath='input_file.xlsx'
    
        appsFacAlumDf = pd.read_excel(inputPath, sheet_name=0)
        appsFacAlumDf.columns = appsFacAlumDf.columns.map(str)

        appsAvailDf = pd.read_excel(inputPath, sheet_name=1)
    
        facAlumAvailDf = pd.read_excel(inputPath, sheet_name=2)
    
        summaryDf = pd.read_excel(inputPath, sheet_name=3)
    
        return appsFacAlumDf, appsAvailDf, facAlumAvailDf, summaryDf
        
    appsFacAlumDf,appsAvailDf,facAlumAvailDf,summaryDf = \
        get_input_file(inputPath)
    
    #Can get the desired counts below if wish:
    #facDesiredCounts  = appsFacAlumDf["AC"].value_counts()
    #Num1DesiredCounts = appsFacAlumDf["1"].value_counts()
    #Num2DesiredCounts = appsFacAlumDf["2"].value_counts()
    
    ################# Get auxiliary input data (.txt) ########################
    ##########################################################################
    ##########################################################################
    
    #Applicant times:
    allTimes = pd.read_table('all_times.txt', 
                             sep='\t', 
                             lineterminator='\n',
                             header=None).replace({'\\r':''}, 
                             regex=True)[0].tolist()

    #Faculty times
    facTimesAll = pd.read_table('all_faculty_time_segments.txt', 
                                sep='\t', 
                                lineterminator='\n',
                                header=None).replace({'\\r':''}, 
                                regex=True)[0].tolist()
    
    #All days available
    allDays = pd.read_table('all_days.txt', 
                            sep='\t', 
                            lineterminator='\n',
                            header=None).replace({'\\r':''}, 
                            regex=True)[0].tolist()
                      
    #All applicants that are available
    apps = pd.read_table('all_applicants.txt', 
                         sep='\t', 
                         lineterminator='\n',
                         header=None).replace({'\\r':''}, 
                         regex=True)[0].tolist()
    
    #Names of all the faculty that are available
    allFaculty = pd.read_table('all_faculty.txt', 
                               sep='\t', 
                               lineterminator='\n',
                               header=None).replace({'\\r':''}, 
                               regex=True)[0].tolist()		
    
    #Names of all alumni that are available
    allAlum = pd.read_table('all_alums.txt', 
                            sep='\t', 
                            lineterminator='\n',
                            header=None).replace({'\\r':''}, 
                            regex=True)[0].tolist()

    #Can combine allFaculty and allAlum DFs if desired:
    #allFacAlum = allFaculty + allAlum

    #Can check num rows if desired:
    #len(allFacAlum) #41
    #len(allTimes) #18
    #len(allDays) #5
    #18*5=90

    ###########################################################################

    #Need to ensure that the input sheet is properly formatted
    #See readme. 
    #For example, we want to be careful about duplicates etc. 
    
    #We can drop unecessary columns or discuss initial input 
    #formatting revisions
    #Below drops everything but names and days
    allAvailFacTimes = facAlumAvailDf.drop(
        facAlumAvailDf.columns[[1,2,3,4,10,11,12,13,14,15]],axis=1)
    allAvailFacTimes = pd.melt(allAvailFacTimes, 
                               id_vars='Your name', 
                               value_vars=allAvailFacTimes.columns[1:],
                               var_name="day",
                               value_name="time")
    allAvailFacTimes.replace('morning', 
                             '''8:00-8:30, 
                             8:30-9:00, 
                             9:00-9:30, 
                             9:30-10:00, 
                             10:00-10:30, 
                             10:30-11:00, 
                             11:00-11:30, 
                             11:30-12:00''', 
                             inplace=True)
    allAvailFacTimes.replace('afternoon',
                             '''12:00-12:30, 
                             12:30-1:00, 
                             1:00-1:30, 
                             1:30-2:00, 
                             2:00-2:30, 
                             2:30-3:00, 
                             3:00-3:30, 
                             3:30-4:00, 
                             4:00-4:30, 
                             4:30-5:00''', 
                             inplace=True)
    allAvailFacTimes.replace('morning, afternoon',
                             '''8:00-8:30, 
                             8:30-9:00, 
                             9:00-9:30, 
                             9:30-10:00, 
                             10:00-10:30, 
                             10:30-11:00, 
                             11:00-11:30, 
                             11:30-12:00, 
                             12:00-12:30, 
                             12:30-1:00, 
                             1:00-1:30, 
                             1:30-2:00, 
                             2:00-2:30, 
                             2:30-3:00, 
                             3:00-3:30, 
                             3:30-4:00, 
                             4:00-4:30, 
                             4:30-5:00''', 
                             inplace=True)

    newTimesFrame = allAvailFacTimes.time.str.split(", ", expand = True)
    allAvailFacTimes = pd.concat([allAvailFacTimes, newTimesFrame], axis=1)
    allAvailFacTimes = allAvailFacTimes.drop(columns='time')
    allAvailFacTimes = pd.melt(allAvailFacTimes, 
                               id_vars=['Your name','day'], 
                               value_vars=allAvailFacTimes.columns[2:],
                               var_name="slot",
                               value_name="time")
    allAvailFacTimes = allAvailFacTimes.drop(columns='slot')
    allAvailFacTimes.loc[:,"applicant"] = "open"
    allAvailFacTimes.columns.values[0] = 'faculty'
    #Drop rows with NA in time:
    allAvailFacTimes.dropna(inplace=True)
    #Drop rows with none in time
    df_filtered = allAvailFacTimes[allAvailFacTimes['time'] != 'none']
    finalDf = allAvailFacTimes.copy()
    finalDf.sort_values(by=['faculty','day','time'],inplace=True,axis=0)
    finalDf.reset_index(inplace=True, drop=True)

    #finalDf['ofaculty'] = pd.Categorical(finalDf.faculty, 
    #                                     ordered=True, 
    #                                     categories=facAlumAvailDf['Your name']
    #                                     )

    #finalDf['oday'] = pd.Categorical(finalDf.day, 
    #                                     ordered=True, 
    #                                     categories=allDays)

    #finalDf['otime'] = pd.Categorical(finalDf.time, 
    #                                     ordered=True, 
    #                                     categories=allTimes)

    #finalDf = finalDf.sort_values(by=['ofaculty','oday'])

    #finalDf.drop('ofaculty', axis=1, inplace=True)
    #finalDf.drop('oday', axis=1, inplace=True)
    #finalDf.drop('otime', axis=1, inplace=True)
    
    #lets make sure there is an entry for every day and time for each faculty
    #if they aren't available we will call it "closed"
    
    etData = np.repeat(allTimes, 
                      len(finalDf['faculty'].unique())*\
                         len(finalDf['day'].unique()))
    nRow = etData.shape[0]
    everyTimeDf = pd.DataFrame({'day':np.repeat(allDays,nRow/len(allDays)),\
                                'time':etData,\
                                'status':np.repeat(['closed'],nRow)})
    everyTimeDf.insert(0, 
                       'faculty', 
                       np.repeat(finalDf['faculty'].unique(),nRow/\
                       len(finalDf['faculty'].unique())),
                       allow_duplicates=True)
    
    finalDf = pd.merge(finalDf, 
                         everyTimeDf,
                         on=['faculty','day','time'],
                         how='right').fillna('closed').drop('status',axis=1)
    finalDf.sort_values(by=['faculty','day','time'],inplace=True,axis=0)
    finalDf.reset_index(inplace=True, drop=True)
    
    ################### Search for faculty 1 ################################
    #########################################################################
    #########################################################################
    #########################################################################
    #########################################################################

    #appsPoints = dict.fromkeys(apps.array,0)
    appsPoints = dict.fromkeys(np.array(apps),0)

    facs = allAvailFacTimes['faculty'].values
    facsAssignNums = dict.fromkeys(facs,0)
    
    facTimes1 = []
    for i in appsPoints.keys():
        #print(i)
        fac1 = appsFacAlumDf.loc[(appsFacAlumDf['Student'] == i),\
                                 ["1"]].iloc[0,0]
        
        #If fac1 isn't listed (is na) or is not included in faculty list
        #select their fac2 as fac1
        if((pd.isna(fac1)) | (facsAssignNums.get(fac1) is None)):
            fac1 = appsFacAlumDf.loc[(appsFacAlumDf['Student'] == i),\
                                     ["2"]].iloc[0,0]
        
        #If fac already has 5, move on
        if(facsAssignNums[fac1]>=5):
            continue
        
        appTimes = appsAvailDf.drop(appsAvailDf\
                                    .columns[[1,7,8,9,10,11,12]],axis=1)\
                                    .loc[(appsAvailDf['Your Name'] == i)]
        
        facTimes = facAlumAvailDf.drop(facAlumAvailDf\
                                       .columns[[1,2,3,4,10,11,12,13,14,15]],
                                       axis=1)\
                                    .loc[(facAlumAvailDf['Your name'] == fac1)]
        
        #Check if first faculty is even available
        if(not facTimes.empty):
            print("faculty is available")
            #now suset days and check first available:
            #print(appTimes)
            #print(facTimes)
            #replace morning and afternoon with list of times above
            #in faculty avail df to do precise matching
            facTimes.replace('morning',
            '8:00-8:30, 8:30-9:00, 9:00-9:30, 9:30-10:00, 10:00-10:30, \
            10:30-11:00, 11:00-11:30, 11:30-12:00', 
            inplace=True)
            
            facTimes.replace('afternoon',
            '12:00-12:30, 12:30-1:00, 1:00-1:30, 1:30-2:00, 2:00-2:30, \
            2:30-3:00, 3:00-3:30, 3:30-4:00, 4:00-4:30, 4:30-5:00', 
            inplace=True)
            
            facTimes.replace('morning, afternoon',
            '8:00-8:30, 8:30-9:00, 9:00-9:30, 9:30-10:00, 10:00-10:30, \
            10:30-11:00, 11:00-11:30, 11:30-12:00, 12:00-12:30, 12:30-1:00, \
            1:00-1:30, 1:30-2:00, 2:00-2:30, 2:30-3:00, 3:00-3:30, 3:30-4:00, \
            4:00-4:30, 4:30-5:00', 
            inplace=True)
            
            #Also can add After5 here ...
            for column in facTimes.columns[1:]:
                #print(facTimes[column])
                #print(appTimes[column])
                if(not facTimes[column].isnull().values.any()):
                    print('fac times present')
                    if(not appTimes.isnull().values.any()):
                        print('app times present')
                        facTimesSplit = facTimes[column].str.split(', ')
                        appTimesSplit = appTimes[column].str.split(', ')
                    
                        facTimesSplit = facTimesSplit.tolist()[0]
                        #Breakup applicant times into half hours:
                        
                        #appTimesSplit = re.sub('8-9am', '8:00-8:30, 8:30-9:00', appTimesSplit)
                        #appTimesSplit = re.sub['9-10am', '9:00-9:30, 9:30-10:00' appTimesSplit]
                        #appTimesSplit = re.sub['10:00-10:30, 10:30-11:00' for "10-11am" in appTimesSplit]
                        #appTimesSplit = re.sub['11:00-11:30, 11:30-12:00' for "11-12pm" in appTimesSplit]
                        #appTimesSplit = re.sub['12:00-12:30, 12:30-1:00' for "12-1pm" in appTimesSplit]
                        #appTimesSplit = re.sub['1:00-1:30, 1:30-2:00' for "1-2pm" in appTimesSplit]
                        #appTimesSplit = re.sub['2:00-2:30, 2:30-3:00' for "2-3pm" in appTimesSplit]
                        #appTimesSplit = re.sub['3:00-3:30, 3:30-4:00' for "3-4pm" in appTimesSplit]
                        #appTimesSplit = re.sub['4:00-4:30, 4:30-5:00' for "4-5pm" in appTimesSplit]
                        
                        appTimesSplit = appTimesSplit.tolist()[0]
                        appTimesSplit = [sub.replace('8-9am', '8:00-8:30, 8:30-9:00') for sub in appTimesSplit]
                        appTimesSplit = [sub.replace('9-10am', '9:00-9:30, 9:30-10:00') for sub in appTimesSplit]
                        appTimesSplit = [sub.replace('10-11am', '10:00-10:30, 10:30-11:00') for sub in appTimesSplit]
                        appTimesSplit = [sub.replace('11-12pm', '11:00-11:30, 11:30-12:00') for sub in appTimesSplit]
                        appTimesSplit = [sub.replace('12-1pm', '12:00-12:30, 12:30-1:00') for sub in appTimesSplit]
                        appTimesSplit = [sub.replace('1-2pm', '1:00-1:30, 1:30-2:00') for sub in appTimesSplit]
                        appTimesSplit = [sub.replace('2-3pm', '2:00-2:30, 2:30-3:00') for sub in appTimesSplit]
                        appTimesSplit = [sub.replace('3-4pm', '3:00-3:30, 3:30-4:00') for sub in appTimesSplit]
                        appTimesSplit = [sub.replace('4-5pm', '4:00-4:30, 4:30-5:00') for sub in appTimesSplit]
                        
                        appTimesSplit = [i.split(', ') for i in appTimesSplit] 
                        #flatten the above list:
                        appTimesSplit = [item for sublist in appTimesSplit for item in sublist]
                        #appTimesSplit.loc[appTimesSplit=="8-9am"] = '8:00-8:30, 8:30-9:00'
                        #appTimesSplit.loc[appTimesSplit=="9-10am"] = '9:00-9:30, 9:30-10:00'
                        #appTimesSplit.loc[appTimesSplit=="10-11am"] = '10:00-10:30, 10:30-11:00'
                        #appTimesSplit.loc[appTimesSplit=="11-12pm"] = '11:00-11:30, 11:30-12:00'
                        #appTimesSplit.loc[appTimesSplit=="12-1pm"] = '12:00-12:30, 12:30-1:00'
                        #appTimesSplit.loc[appTimesSplit=="1-2pm"] = '1:00-1:30, 1:30-2:00'
                        #appTimesSplit.loc[appTimesSplit=="2-3pm"] = '2:00-2:30, 2:30-3:00'
                        #appTimesSplit.loc[appTimesSplit=="3-4pm"] = '3:00-3:30, 3:30-4:00'
                        #appTimesSplit.loc[appTimesSplit=="4-5pm"] = '4:00-4:30, 4:30-5:00'
                    
                        #look for intersection and take first instance of intersection
                        olap = list(set(facTimesSplit) & set(appTimesSplit))
                        #print(olap)
                        if(bool(olap)==True):
                            print('overlap between facTimes and appTimes')
                            #Check if slot is already assigned in final sheet,
                            #if not, Assign the first available overlap to the final sheet
                            for check in range(len(olap)):
                                #exit()
    ############## major if case here ###########################################                            
                                #check present time slot:
                                checkPresent = (finalDf[(finalDf['faculty']==fac1) & (finalDf['time']==olap[check]) & (finalDf['day']==column)]['applicant'].to_string(index=False) == 'open') & (appsPoints[i] < 2)
                                #get row index of current time slot
                                currentIndex = finalDf.index[(finalDf['faculty']==fac1) & (finalDf['time']==olap[check]) & (finalDf['day']==column)].values[0]
                                #check 1/2 hour after 8 and 1/2 hour before 5:
                                checkBefore = (not ((finalDf.iloc[currentIndex]['applicant'] == '8:00-8:30') & (finalDf.iloc[currentIndex+1]['applicant'] != 'open')))
                                checkAfter = (not ((finalDf.iloc[currentIndex]['applicant'] == '4:30-5:00') & (finalDf.iloc[currentIndex-1]['applicant'] != 'open')))
                                checkInBetween = ((finalDf.iloc[currentIndex]['applicant'] != '8:00-8:30') & (finalDf.iloc[currentIndex]['applicant'] != '4:30-5:00') & (finalDf.iloc[currentIndex-1]['applicant'] == 'open') & (finalDf.iloc[currentIndex+1]['applicant'] == 'open'))
                                if((checkPresent) & (checkBefore) & (checkAfter) & (checkInBetween)):
                                    print('assigning')
                                    #exit()
                                    #print(olap[check])
                                    finalDf.loc[(finalDf['faculty']==fac1) &\
                                            (finalDf['time']==olap[check]) &\
                                            (finalDf['day']==column),
                                            ['applicant']] = i
                                    appsPoints[i]+=1
                                    facsAssignNums[fac1]+=1
                                    selectDf = finalDf.loc[(finalDf['faculty']==fac1) &\
                                            (finalDf['time']==olap[check]) &\
                                            (finalDf['day']==column) &\
                                            (finalDf['applicant']==i)]
                                    print(selectDf['applicant'].to_string(index=False) + ',\t' +\
                                            selectDf['faculty'].to_string(index=False) + ',\t' +\
                                            selectDf['day'].to_string(index=False) + ',\t' +\
                                            selectDf['time'].to_string(index=False))
                                    vals = selectDf.astype(str).values.flatten().tolist()
                                    stringOut = ';\t'.join(map(str,vals))
                                    facTimes1.append(stringOut)
                                    break
                            else:
                                continue
                            break
                             
                                #for check in range(len(olap)):
                                #    if(any(finalDf['faculty'].str.match(fac1)) & any((finalDf['time'].str.match(olap[check]))) & any((finalDf['day'].str.match(column))) & any((finalDf['applicant']==i))):
                                #        #print('already assigned')
                                #        next
                                #    else:
                                #        print(olap[check])
                                #        finalDf.loc[(finalDf['faculty']==fac1) &\
                                #            (finalDf['time']==olap[check]) &\
                                #            (finalDf['day']==column),
                                #            ['applicant']] = i
                                #        appsPoints[i]+=1
                                #        break

    #Output fac1 matches file:
    if os.path.exists("fac1.txt"):
        os.remove("fac1.txt")

    with open('fac1.txt', 'w') as f1:
        f1.write('\n'.join([str(x) for x in facTimes1]))
    f1.close()    
    
########################### Search for faculty 2 ###########################
############################################################################
############################################################################
############################################################################
############################################################################

#print("\n\n\n\n\n fac 2 \n\n\n\n\n ")

    facTimes2 = []
    for i in appsPoints.keys():
        #print(i)
        #selected = [appsFacAlumDf['Student'==i,'1'],appsFacAlumDf['Student'==i,"2"]]
        #appsList.append(appsFacAlumDf.loc[(appsFacAlumDf['Student'] == i),["1","2"]])
        #Find first faculty:
        fac2 = appsFacAlumDf.loc[(appsFacAlumDf['Student'] == i),["2"]].iloc[0,0]
    
        #If fac1 isn't listed (is na) or is not included in faculty list
        #select their fac2 as fac3 or just skip
        if((pd.isna(fac2)) | (facsAssignNums.get(fac2) is None)):
            #fac2 = appsFacAlumDf.loc[(appsFacAlumDf['Student'] == i),["3"]].iloc[0,0]
            continue
    
        #If fac already has 5, move on
        if(facsAssignNums[fac2]>=5):
            continue
    
        appTimes = appsAvailDf.drop(appsAvailDf.columns[[1,7,8,9,10,11,12]],axis=1).loc[(appsAvailDf['Your Name'] == i)]
        facTimes = facAlumAvailDf.drop(facAlumAvailDf.columns[[1,2,3,4,10,11,12,13,14,15]],axis=1).loc[(facAlumAvailDf['Your name'] == fac2)]
        #Check if first faculty is even available
        if(not facTimes.empty):
            #now suset days and check first available:
            #print(appTimes)
            #print(facTimes)
            #replace morning and afternoon with list of times above
            #in faculty avail df to do precise matching
            facTimes.replace('morning', '8:00-8:30, 8:30-9:00, 9:00-9:30, 9:30-10:00, 10:00-10:30, 10:30-11:00, 11:00-11:30, 11:30-12:00', inplace=True)
            facTimes.replace('afternoon','12:00-12:30, 12:30-1:00, 1:00-1:30, 1:30-2:00, 2:00-2:30, 2:30-3:00, 3:00-3:30, 3:30-4:00, 4:00-4:30, 4:30-5:00', inplace=True)
            facTimes.replace('morning, afternoon','8:00-8:30, 8:30-9:00, 9:00-9:30, 9:30-10:00, 10:00-10:30, 10:30-11:00, 11:00-11:30, 11:30-12:00, 12:00-12:30, 12:30-1:00, 1:00-1:30, 1:30-2:00, 2:00-2:30, 2:30-3:00, 3:00-3:30, 3:30-4:00, 4:00-4:30, 4:30-5:00', inplace=True)
            #Also can add After5 here ...
            for column in facTimes.columns[1:]:
                #print(facTimes[column])
                #print(appTimes[column])
                if(not facTimes[column].isnull().values.any()):
                    #check if app time is NA
                    if(not appTimes.isnull().values.any()):
                        facTimesSplit = facTimes[column].str.split(', ')
                        appTimesSplit = appTimes[column].str.split(', ')
                
                        facTimesSplit = facTimesSplit.tolist()[0]
                        #Breakup applicant times into half hours:
                    
                        appTimesSplit = appTimesSplit.tolist()[0]
                        appTimesSplit = [sub.replace('8-9am', '8:00-8:30, 8:30-9:00') for sub in appTimesSplit]
                        appTimesSplit = [sub.replace('9-10am', '9:00-9:30, 9:30-10:00') for sub in appTimesSplit]
                        appTimesSplit = [sub.replace('10-11am', '10:00-10:30, 10:30-11:00') for sub in appTimesSplit]
                        appTimesSplit = [sub.replace('11-12pm', '11:00-11:30, 11:30-12:00') for sub in appTimesSplit]
                        appTimesSplit = [sub.replace('12-1pm', '12:00-12:30, 12:30-1:00') for sub in appTimesSplit]
                        appTimesSplit = [sub.replace('1-2pm', '1:00-1:30, 1:30-2:00') for sub in appTimesSplit]
                        appTimesSplit = [sub.replace('2-3pm', '2:00-2:30, 2:30-3:00') for sub in appTimesSplit]
                        appTimesSplit = [sub.replace('3-4pm', '3:00-3:30, 3:30-4:00') for sub in appTimesSplit]
                        appTimesSplit = [sub.replace('4-5pm', '4:00-4:30, 4:30-5:00') for sub in appTimesSplit]
                    
                        appTimesSplit = [i.split(', ') for i in appTimesSplit] 
                        #flatten the above list:
                        appTimesSplit = [item for sublist in appTimesSplit for item in sublist]
                
                        #look for intersection and take first instance of intersection
                        olap = list(set(facTimesSplit) & set(appTimesSplit))
                        #print(olap)
                        if(bool(olap)==True):
                            #Check if slot is already assigned in final sheet,
                            #if not, Assign the first available overlap to the final sheet
                            for check in range(len(olap)):
                                #exit()
############## major if case here ###########################################                            

                                #check present time slot:
                                checkPresent = (finalDf[(finalDf['faculty']==fac2) & (finalDf['time']==olap[check]) & (finalDf['day']==column)]['applicant'].to_string(index=False) == 'open') & (appsPoints[i] < 2)
                                #get row index of current time slot
                                currentIndex = finalDf.index[(finalDf['faculty']==fac2) & (finalDf['time']==olap[check]) & (finalDf['day']==column)][0]
                                #check 1/2 hour after 8 and 1/2 hour before 5:
                                checkBefore = (not ((finalDf.iloc[currentIndex]['applicant'] == '8:00-8:30') & \
                                    (finalDf.iloc[currentIndex+1]['applicant'] != 'open')))
                                checkAfter = (not ((finalDf.iloc[currentIndex]['applicant'] == '4:30-5:00') & \
                                    (finalDf.iloc[currentIndex-1]['applicant'] != 'open')))
                                checkInBetween = ((finalDf.iloc[currentIndex]['applicant'] != '8:00-8:30') & \
                                              (finalDf.iloc[currentIndex]['applicant'] != '4:30-5:00') & \
                                              (finalDf.iloc[currentIndex-1]['applicant'] == 'open') & \
                                              (finalDf.iloc[currentIndex+1]['applicant'] == 'open'))
                                
                                if((checkPresent) & (checkBefore) & (checkAfter) & (checkInBetween)):
                                    print('assigning')
                                    #exit()
                                    #print(olap[check])
                                    finalDf.loc[(finalDf['faculty']==fac2) &\
                                        (finalDf['time']==olap[check]) &\
                                        (finalDf['day']==column),
                                        ['applicant']] = i
                                    appsPoints[i]+=1
                                    facsAssignNums[fac2]+=1
                                    selectDf = finalDf.loc[(finalDf['faculty']==fac2) &\
                                        (finalDf['time']==olap[check]) &\
                                        (finalDf['day']==column) &\
                                        (finalDf['applicant']==i)]
                                    print(selectDf['applicant'].to_string(index=False) + ',\t' +\
                                        selectDf['faculty'].to_string(index=False) + ',\t' +\
                                        selectDf['day'].to_string(index=False) + ',\t' +\
                                        selectDf['time'].to_string(index=False))
                                    vals = selectDf.astype(str).values.flatten().tolist()
                                    stringOut = ';\t'.join(map(str,vals))
                                    facTimes2.append(stringOut)
                                    break
                        else:
                            continue
                        break

    #Output fac2 matches file:
        if os.path.exists("fac2.txt"):
            os.remove("fac2.txt")
    
    with open('fac2.txt', 'w') as f2:
        f2.write('\n'.join([str(x) for x in facTimes2]))
        f2.close()  
    
######### write out applicant points file and fac #s files ############
#######################################################################
#######################################################################
#######################################################################
#######################################################################

    dfp = pd.DataFrame(data=appsPoints, index=[0])
    #Transpose:
    dfp = (dfp.T)
    dfpIndex = dfp.index
    dfp.insert(0,"applicant",dfpIndex)
    dfp.columns.values[1] = 'points'
    dfp.to_excel('apps_points.xlsx',index=False)

    dfa = pd.DataFrame(data=facsAssignNums, index=[0])
    #Transpose:
    dfa = (dfa.T)
    dfaIndex = dfa.index
    dfa.insert(0,"faculty",dfaIndex)
    dfa.columns.values[1] = 'assigns'
    dfa.to_excel('fac_num_assigned.xlsx',index=False)    
    
    ################ Assigning remaining slots randomly ###################
    #######################################################################
    #######################################################################
    #######################################################################
    #######################################################################

    #remember to prioritize faculty over alumni
    facSubs = finalDf[~finalDf['faculty'].isin(allAlum)]['faculty'].unique()
    alumSubs = finalDf[finalDf['faculty'].isin(allAlum)]['faculty'].unique()
    allSubs = facSubs.tolist()+alumSubs.tolist()
    allJobs = (np.repeat('fac',len(facSubs)).tolist()+np.repeat('alum',len(alumSubs)).tolist())

    alumFacCheck = pd.DataFrame()
    alumFacCheck['faculty']=allSubs
    alumFacCheck['job']=allJobs
    alumFacCheck['status']=np.repeat(0,alumFacCheck.shape[0])
                                 
    assignRandom = True
    if(assignRandom==True):
        #reindex row #s:
            finalDf.reset_index(inplace=True, drop=True)                            
            #Now loop through and fill in missing
            for i in appsPoints.keys():
                print("assigning random to " + i)
                #print(len(finalDf[finalDf['applicant']==i]))
                #if no 2 points, choose random faculty/alum for the candidate
                while(appsPoints[i]<2):
                    checkPresent = False
                    checkBefore = False
                    checkAfter = False
                    checkInBetween = False
                    #pic random row in assignment df:
                    randomRow = finalDf.sample()
                    #print(randomRow['faculty'])
                    #check if faculty:
                    if(((alumFacCheck.loc[alumFacCheck['faculty']==randomRow.faculty.values[0],'job'].values[0]=='fac') &\
                                (alumFacCheck.loc[alumFacCheck['faculty']==randomRow.faculty.values[0],'status'].values[0]==0)) |\
                               ((all(alumFacCheck.loc[alumFacCheck['job']=='fac'].status==1) & \
                                 (alumFacCheck.loc[alumFacCheck['faculty']==randomRow.faculty.values[0],'job'].values[0]=='alum')))):
                                
                        #check present time slot:
                        checkPresent = (randomRow['applicant'].to_string(index=False)=='open')
                                    
                        #get row index of current time slot
                        currentIndex = randomRow.index[0]
            
                        #check fac #s alreads assigned:
                        selFac = finalDf.iloc[currentIndex]['faculty']
                
                        #assign the faculty over alum:
                        alumFacCheck.loc[alumFacCheck['faculty']==selFac,'status'] = 1
                
                        if(facsAssignNums[selFac]>=5):
                            continue
            
                        #check if index is first or last
                        if(currentIndex==len(finalDf.index)-1):
                            checkBefore = (finalDf.iloc[currentIndex-1]['applicant'] != 'open')
                            checkAfter = True
                            checkInBetween = True
                        elif(currentIndex==0):
                            checkAfter = (finalDf.iloc[currentIndex+1]['applicant'] != 'open')
                            checkBefore = True
                            checkInBetween = True
                        else:
                            #check 1/2 hour after 8 and 1/2 hour before 5:
                            checkBefore = (not ((finalDf.iloc[currentIndex]['applicant'] == '8:00-8:30') & (finalDf.iloc[currentIndex+1]['applicant'] != 'open')))
                            checkAfter = (not ((finalDf.iloc[currentIndex]['applicant'] == '4:30-5:00') & (finalDf.iloc[currentIndex-1]['applicant'] != 'open')))
                            checkInBetween = ((finalDf.iloc[currentIndex]['applicant'] != '8:00-8:30') & (finalDf.iloc[currentIndex]['applicant'] != '4:30-5:00') & (finalDf.iloc[currentIndex-1]['applicant'] == 'open') & (finalDf.iloc[currentIndex+1]['applicant'] == 'open'))
                            #print(checkPresent + checkBefore + checkAfter + checkInBetween)
                        if((checkPresent) & (checkBefore) & (checkAfter) & (checkInBetween)):
                            #print('assigning random')
                            finalDf.iloc[currentIndex]["applicant"] = i 
                            appsPoints[i]+=1
                            facsAssignNums[selFac]+=1
                            alumFacCheck.loc[alumFacCheck['faculty']==selFac,'status'] = 1
                    
                            #Reset faculty over alum priority check:
                            alumFacCheck['faculty']=allSubs
                            alumFacCheck['job']=allJobs
                            alumFacCheck['status']=np.repeat(0,alumFacCheck.shape[0])
                    
                            #print(appsPoints[i])
                            #print(appsPoints)

    print(alumFacCheck)   
    
################### Final formatting and output #######################
#######################################################################
#######################################################################
#######################################################################
#######################################################################

    finalDups = finalDf.duplicated()

    finalDf = finalDf.drop_duplicates(finalDf.columns).reset_index(drop=True)
    
    finalDfM  = finalDf.loc[finalDf["day"]==allDays[0]]
    finalDfM = finalDfM.drop_duplicates(finalDfM.columns).reset_index(drop=True)

    finalDfT  = finalDf.loc[finalDf["day"]==allDays[1]] 
    finalDfT = finalDfT.drop_duplicates(finalDfT.columns).reset_index(drop=True)

    finalDfW  = finalDf.loc[finalDf["day"]==allDays[2]]
    finalDfW = finalDfW.drop_duplicates(finalDfW.columns).reset_index(drop=True)

    finalDfTh = finalDf.loc[finalDf["day"]==allDays[3]] 
    finalDfTh = finalDfTh.drop_duplicates(finalDfTh.columns).reset_index(drop=True)

    finalDfF  = finalDf.loc[finalDf["day"]==allDays[4]] 
    finalDfF = finalDfF.drop_duplicates(finalDfF.columns).reset_index(drop=True)

    len(finalDfM["faculty"].unique()) #27
    len(finalDfT["faculty"].unique()) #23
    len(finalDfW["faculty"].unique()) #24
    len(finalDfTh["faculty"].unique()) #20
    len(finalDfF["faculty"].unique()) #24
    
    #Need to remove blank time spot for a couple places:
    finalDfW.drop([4], axis=0, inplace=True)
    finalDfF.drop([6], axis=0, inplace=True)

    mTimes = finalDfM['time']
    finalDfMCast  = (finalDfM.pivot(index=['time'], 
                                         columns='faculty',
                                         values='applicant').reset_index(drop=True))
    finalDfMCast.insert(0,'time',mTimes.unique())
    tTimes = finalDfT['time']
    finalDfTCast  = (finalDfT.pivot(index=['time','day'], 
                                         columns='faculty',
                                         values='applicant').reset_index(drop=True))
    finalDfTCast.insert(0,'time',tTimes.unique())
    wTimes = finalDfW['time']
    finalDfWCast  = (finalDfW.pivot(index=['time','day'], 
                                         columns='faculty',
                                         values='applicant').reset_index(drop=True))
    finalDfWCast.insert(0,'time',wTimes.unique())
    thTimes = finalDfTh['time']
    finalDfThCast = (finalDfTh.pivot(index=['time','day'], 
                                         columns='faculty',
                                         values='applicant').reset_index(drop=True))
    finalDfThCast.insert(0,'time',thTimes.unique())
    fTimes = finalDfF['time']
    finalDfFCast  = (finalDfF.pivot(index=['time','day'], 
                                         columns='faculty',
                                         values='applicant').reset_index(drop=True))
    finalDfFCast.insert(0,'time',fTimes.unique())

    len(finalDfMCast)  #18 rows
    len(finalDfTCast)  #18 rows
    len(finalDfWCast)  #19 rows
    len(finalDfThCast) #18 rows
    len(finalDfFCast)  #19 rows

    #finalDfWCast= finalDfWCast.drop([18],axis=0)
    #finalDfFCast= finalDfFCast.drop(18,axis=0)
    
    finalDfMCast.fillna("", inplace=True)
    finalDfTCast.fillna("", inplace=True)
    finalDfWCast.fillna("", inplace=True)
    finalDfThCast.fillna("", inplace=True)
    finalDfFCast.fillna("", inplace=True)

    finalDfMCast.replace('open', "", inplace=True)
    finalDfTCast.replace('open', "", inplace=True)
    finalDfWCast.replace('open', "", inplace=True)
    finalDfThCast.replace('open', "", inplace=True)
    finalDfFCast.replace('open', "", inplace=True)

    #unique(finalDfM["time"])
    #unique(finalDfT["time"])
    #unique(finalDfW["time"])
    #unique(finalDfTh["time"])
    #unique(finalDfF["time"])

    #finalDfMCast.insert(0,finalDfM inplace=TRUE

    #Order by time column:
    finalDfMCast['otime'] = pd.Categorical(finalDfMCast.time, 
                                       ordered=True, 
                                       categories=allTimes)
    finalDfMCast = finalDfMCast.sort_values('otime')
    finalDfMCast.drop('otime',axis=1,inplace=True)

    finalDfTCast['otime'] = pd.Categorical(finalDfTCast.time, 
                                       ordered=True, 
                                       categories=allTimes)
    finalDfTCast = finalDfTCast.sort_values('otime')
    finalDfTCast.drop('otime',axis=1,inplace=True)

    finalDfWCast['otime'] = pd.Categorical(finalDfWCast.time, 
                                       ordered=True, 
                                       categories=allTimes)
    finalDfWCast = finalDfWCast.sort_values('otime')
    finalDfWCast.drop('otime',axis=1,inplace=True)

    finalDfThCast['otime'] = pd.Categorical(finalDfThCast.time, 
                                       ordered=True, 
                                       categories=allTimes)
    finalDfThCast = finalDfThCast.sort_values('otime')
    finalDfThCast.drop('otime',axis=1,inplace=True)

    finalDfFCast['otime'] = pd.Categorical(finalDfFCast.time, 
                                       ordered=True, 
                                       categories=allTimes)
    finalDfFCast = finalDfFCast.sort_values('otime')
    finalDfFCast.drop('otime',axis=1,inplace=True)

    #Now add the day info to the data frames as second column
    #and row bind the casted table and output to file
    #del finalDfMCast
    #del finalDfTCast
    #del finalDfWCast
    #del finalDfThCast
    #del finalDfFCast
    
    dfList=[finalDfMCast,
            finalDfTCast,
            finalDfWCast,
            finalDfThCast,
            finalDfFCast]
    
    #dfList=[finalDfMCast.insert(loc=1, column="day", value = allDays[0]),
    #        finalDfTCast.insert(loc=1, column="day", value = allDays[1]),
    #        finalDfWCast.insert(loc=1, column="day", value = allDays[2]),
    #        finalDfThCast.insert(loc=1, column="day", value = allDays[3]),
    #        finalDfFCast.insert(loc=1, column="day", value = allDays[4])]
    
    finalOutputFrame = pd.concat(dfList, axis=0)
    #finalOutputFrame = pd.concat(dfList, axis=1)
    #outputDataTable = dt.Frame.from_pandas(finalOutputFrame)
    
################# paste points to the applicant ##################
################# paste job and assign num to faculty/alum #######

    #Write to excel:
    writer = pd.ExcelWriter("final_match_output.xlsx",engine='xlsxwriter') 
    startRow = 2
    headerRow = 0
    workbook=writer.book
    count = 0
    for dataFrame in dfList:
        dataFrame.style.set_properties(**{'background-color': 'blue'}, subset=['A'])
        dataFrame.to_excel(writer, 
                           sheet_name='output', 
                           startrow=startRow, 
                           startcol=0,
                           index=False)
        bold = workbook.add_format({'bold': True})
        worksheet = writer.sheets['output']
        worksheet.write(headerRow, 0, allDays[count],bold)
        #worksheet.conditional_format('B2:B8', {'type': '3_color_scale'})
        startRow+=22
        headerRow+=22
        count+=1
    writer.close()   
    
######################### Server Action ##################################
##########################################################################
##########################################################################

    mod_counter = reactive.value(0)

    query_output_server("initial_query", con=con, remove_id="initial_query")
    
    @reactive.effect
    @reactive.event(input.add_query)
    def _():
        counter = mod_counter.get() + 1
        mod_counter.set(counter)
        id = "query_" + str(counter)
        ui.insert_ui(
            selector="#module_container",
            where="afterBegin",
            ui=query_output_ui(id, remove_id=id),
        )
        query_output_server(id, con=con, remove_id=id)
    
    @reactive.effect
    @reactive.event(input.show_meta)
    def _():
        counter = mod_counter.get() + 1
        mod_counter.set(counter)
        id = "query_" + str(counter)
        ui.insert_ui(
            selector="#module_container",
            where="afterBegin",
            ui=query_output_ui(
                id, qry="SELECT * from information_schema.columns", remove_id=id
            ),
        )
        query_output_server(id, con=con, remove_id=id)

    @render.text
    #@reactive.event(input.radio_group)
    def radio_value():
        return f"You selected: {input.radio_group()}"
    
    @render.text
    #@reactive.event(input.radio_group)
    def radio_frame():
        match input.radio_group():
            case 'Show applicants':
                return '\n'.join(str(x) for x in apps)#.replace(r'\n', '\n')
            case  'Show faculty':
                return '\n'.join(str(x) for x in allFaculty)
            case  'Show alumni':
                return '\n'.join(str(x) for x in allAlum)
            case  'Show times':
                return '\n'.join(str(x) for x in allTimes)
            case  'Show faculty times':
                return '\n'.join(str(x) for x in facTimesAll)
            case  'Show days':
                return '\n'.join(str(x) for x in allDays)
        #return f"You selected: {input.radio_group()}"
        #appTimesAll
        #facTimesAll
        #apps
        #allFaculty
        #allAlum
        #allTimes    
        #allDays
    
    #df: reactive.value[pd.DataFrame] = reactive.value(
    #    #sns.load_dataset("anagrams").iloc[:, 1:]
    #    finalOutputFrame
    #)
    
    #height = 350
    #width = "fit-content"
    @render.data_frame
    @reactive.event(input.run_matches)
    def table():
        #return render.DataTable(
        #    #df(),
        #    finalOutputFrame,
        #    width=width,
        #    height=height,
        #    filters=input.filters(),
        #    editable=input.editable(),
        #    selection_mode=input.selection_mode(),
        return finalOutputFrame
        #)  
        
    @render.download()
    def downloadData():
        """
        This is the simplest case. The implementation simply returns the name of a file.
        Note that the function name (`download1`) determines which download_button()
        corresponds to this function.
        """

        path = os.path.join(os.path.dirname(__file__), "final_match_output.xlsx")
        return path
     
app = App(app_ui, server)

######### For running the App in local browser ###########################
app.run()
##########################################################################
