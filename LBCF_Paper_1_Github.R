#load libraries
library(Rcpp)
library(bcf)
library(dbarts)
library(missRanger)
library(naniar)
library(dplyr)
library(tidyr)
library(plyr)
library(mltools)
library(data.table)
library(bartMachine)


#load data
load("~/your_directory/HSLS_STUDENT_DATA.rdata")
df<-hsls_17_student_pets_sr_v1_0
rm(hsls_17_student_pets_sr_v1_0)


###################################
###################################
#Wave 1 variables
#Change X1TXMTH1 to X1TXMTH2 etc 
#for other plausible values
###################################
###################################


mydf1<-df[,c("X1SEX",
             "X1RACE",
             "X1DUALLANG",
             "X1STDOB",
             "X1TXMTH1",
             "X1PARPATTERN",
             "X1PAREDU",
             "X1HHNUMBER",
             "X1SES",
             "X1MTHID",
             "X1MTHUTI",
             "X1MTHEFF",
             "X1MTHINT",
             "X1SCHOOLBEL",
             "X1SCHOOLENG",
             "X1STU30OCC_STEM1",
             "X1STUEDEXPCT",
             "X1CONTROL",
             "X1LOCALE",
             "X1REGION",
             "X1SCHOOLCLI",
             "S1MCLUB",
             "S1MCOMPETE",
             "S1MCAMP",
             "S1MTUTOR",
             "S1MFALL09",
             "S1HRMHOMEWK",
             "S1HRWORK",
             "X2TXMTH1",
             "STU_ID",
             "S1M8",
             "S1M8GRADE",
             "S1S8GRADE",
             "W1STUDENT")]

mydf1$STU_ID<-as.numeric(as.character(mydf1$STU_ID))

mydf1<-mydf1[order(mydf1$STU_ID),]

mydf1$X1TXMTH1<-ifelse(mydf1$X1TXMTH1<(-6), NA, mydf1$X1TXMTH1)

ach1_not_valid<-is.na(mydf1$X1TXMTH1)

mydf1<-mydf1[!ach1_not_valid,]

mydf1$X1SEX<-ifelse(mydf1$X1SEX=="Male", 1, 0)

mydf1$X1DUALLANG<-ifelse(mydf1$X1DUALLANG=="Missing", NA, mydf1$X1DUALLANG)

age_fun<-function(agenum)
{
  if(agenum=="-9")
  {
    return(NA)
  }
  yearnum<-as.numeric(as.character(substr(agenum, 1, 4)))
  monthnum<-as.numeric(as.character(substr(agenum, 5, 6)))
  months_old<-(2009-yearnum)*12-monthnum
  return(months_old)
}

mydf1$X1STDOB<-as.numeric(lapply(mydf1$X1STDOB, age_fun))

from<-c("Less than high school", "High school diploma or GED", 
        "Associate's degree", "Bachelor's degree", "Master's degree",                         
        "Ph.D/M.D/Law/other high lvl prof degree", "Unit non-response")

to<-c(1, 2, 3, 4, 5, 6, 0)

mymap<-function(x)
{
  x<-mapvalues(x,
               from=from,
               to=to)
  
  x<-as.numeric(as.character(x))
  
  return(x)
}

#Best way to habdle this?
mydf1<-mutate_at(mydf1, c("X1PAREDU"), mymap)

mydf1$parent_missing<-ifelse(mydf1$X1PAREDU==0, 1, 0)

hh_fun<-function(hhnum)
{
  if(hhnum=="Unit non-response")
  {
    return(0)
  }
  if(hhnum=="Missing")
  {
    return(NA)
  }
  num<-as.numeric(as.character(substr(hhnum, 1, 2)))
  return(num)
}

mydf1$X1HHNUMBER<-as.numeric(lapply(mydf1$X1HHNUMBER, hh_fun))

mydf1$X1MTHID<-ifelse(mydf1$X1MTHID==-9, NA, mydf1$X1MTHID)
mydf1$X1MTHUTI<-ifelse(mydf1$X1MTHUTI==-9, NA, mydf1$X1MTHUTI)
mydf1$X1MTHEFF<-ifelse(mydf1$X1MTHEFF==-9, NA, mydf1$X1MTHEFF)
mydf1$X1MTHINT<-ifelse(mydf1$X1MTHIN==-9, NA, mydf1$X1MTHINT)
mydf1$X1SCHOOLBEL<-ifelse(mydf1$X1SCHOOLBEL==-9, NA, mydf1$X1SCHOOLBEL)
mydf1$X1SCHOOLENG<-ifelse(mydf1$X1SCHOOLENG==-9, NA, mydf1$X1SCHOOLENG)

mydf1$X1STU30OCC_STEM1<-ifelse(mydf1$X1STU30OCC_STEM1=="Not a STEM occupation", 0, 
                               ifelse(mydf1$X1STU30OCC_STEM1=="Missing", NA, 1))

from<-c("Less than high school", "High school diploma or GED",             
        "Start an Associate's degree", "Complete an Associate's degree",         
        "Start a Bachelor's degree", "Complete a Bachelor's degree",           
        "Start a Master's degree", "Complete a Master's degree",             
        "Start Ph.D/M.D/Law/other prof degree", "Complete Ph.D/M.D/Law/other prof degree",
        "Don't know")

to<-c(1, 2, 3, 3, 4, 4, 5, 5, 6, 6, 0)

mydf1<-mutate_at(mydf1, c("X1STUEDEXPCT"), mymap)

mydf1$X1CONTROL<-ifelse(mydf1$X1CONTROL=="Public", 1, 0)

mydf1$no_principal<-ifelse(mydf1$X1SCHOOLCLI==-8, 1, 0)

mydf1$X1SCHOOLCLI<-ifelse(mydf1$X1SCHOOLCLI==-9, NA, mydf1$X1SCHOOLCLI)

mydf1$S1MCLUB<-ifelse(mydf1$S1MCLUB=="Yes", 1, 
                      ifelse(mydf1$S1MCLUB=="Missing", NA, 0))

mydf1$S1MCOMPETE<-ifelse(mydf1$S1MCOMPETE=="Yes", 1, 
                         ifelse(mydf1$S1MCOMPETE=="Missing", NA, 0))

mydf1$S1MCAMP<-ifelse(mydf1$S1MCAMP=="Yes", 1, 
                      ifelse(mydf1$S1MCAMP=="Missing", NA, 0))

mydf1$S1MTUTOR<-ifelse(mydf1$S1MTUTOR=="Yes", 1, 
                       ifelse(mydf1$S1MTUTOR=="Missing", NA, 0))

mydf1$S1MFALL09<-ifelse(mydf1$S1MFALL09=="Yes", 1, 
                        ifelse(mydf1$S1MFALL09=="Missing", NA, 0))


from<-c("Less than 1 hour", "1 to 2 hours", "2 to 3 hours",           
        "3 to 4 hours", "4 to 5 hours", "5 or more hours",        
        "Item legitimate skip/NA", "Missing")

to<-c(1, 2, 3, 4, 5, 6, 0, NA)

mydf1<-mutate_at(mydf1, c("S1HRMHOMEWK"), mymap)

mydf1<-mutate_at(mydf1, c("S1HRWORK"), mymap)

from<-c("Math 8", "Advanced or Honors Math 8", "Pre-algebra",
        "Algebra I including IA and IB", "Algebra II or Trigonometry", "Geometry", 
        "Integrated Math", "Other math course", "Item legitimate skip/NA", 
        "Unit non-response", "Missing")

to<-c(1, 2, 1, 2, 3, 3, 1, 1, 0, 0, NA)

mydf1<-mutate_at(mydf1, c("S1M8"), mymap)

from<-c("A", "B", "C", "D", "Below D", "Class was not graded",
        "Item legitimate skip/NA", "Unit non-response", "Missing")

to<-c(5,4,3,2,1,0,0,0,NA)

mydf1<-mutate_at(mydf1, c("S1M8GRADE"), mymap)

mydf1<-mutate_at(mydf1, c("S1S8GRADE"), mymap)



###################################
###################################
#Wave 2 variables
###################################
###################################


mydf2<-df[,c("X2SEX",
             "X2RACE",
             "X2DUALLANG",
             "X2STDOB",
             "X2TXMTH1",
             "X2PARPATTERN",
             "X2PAREDU",
             "X2HHNUMBER",
             "X2SES",
             "X2MTHID",
             "X2MTHINT_R",
             "X2MTHUTI",
             "X2MTHEFF",
             "X2STU30OCC_STEM1",
             "X2STUEDEXPCT",
             "X2CONTROL",
             "X2LOCALE",
             "X2REGION",
             "X2SCHOOLCLI",
             "S2MCLUB",
             "S2MCOMPETE",
             "S2MHOMEWRK",
             "S2HSJOBHR",
             "S2MTUTORED",
             "X2EVERDROP",
             "X2BEHAVEIN",
             "X2MEFFORT",
             "X2PROBLEM",
             "S2APMATH",
             "S2APSCIENCE",
             "S2MSPR12",
             "S2HSJOBEVER",
             "X1TXMTH1",
             "STU_ID",
             "S2HIMATH12",
             "X2NUMHS",
             "W2STUDENT")]

mydf2$STU_ID<-as.numeric(as.character(df$STU_ID))

mydf2<-mydf2[order(mydf2$STU_ID),]

mydf2$X2TXMTH1<-ifelse(mydf2$X2TXMTH1<(-6), NA, mydf2$X2TXMTH1)

ach2_not_valid<-is.na(mydf2$X2TXMTH1)

mydf2<-mydf2[!ach1_not_valid & !ach2_not_valid,]

mydf2$X2SEX<-ifelse(mydf2$X2SEX=="Male", 1, 0)

mydf2$X2DUALLANG<-ifelse(mydf2$X2DUALLANG=="Missing", NA, mydf2$X2DUALLANG)

age_fun<-function(agenum)
{
  if(agenum=="-9")
  {
    return(NA)
  }
  yearnum<-as.numeric(as.character(substr(agenum, 1, 4)))
  monthnum<-as.numeric(as.character(substr(agenum, 5, 6)))
  months_old<-(2009-yearnum)*12-monthnum
  return(months_old)
}

mydf2$X2STDOB<-as.numeric(lapply(mydf2$X2STDOB, age_fun))

from<-c("Less than high school",                                          
        "High school diploma or GED or alterntive HS credential",         
        "Certificate/diploma from school providing occupational training",
        "Associate's degree",                                             
        "Bachelor's degree",                                              
        "Master's degree",                                                
        "Ph.D/M.D/Law/other high lvl prof degree")

to<-c(1, 2, 2, 3, 4, 5, 6)

mymap<-function(x)
{
  x<-mapvalues(x,
               from=from,
               to=to)
  
  x<-as.numeric(as.character(x))
  
  return(x)
}

mydf2<-mutate_at(mydf2, c("X2PAREDU"), mymap)

hh_fun<-function(hhnum)
{
  if(hhnum=="Unit non-response")
  {
    return(0)
  }
  if(hhnum=="Missing")
  {
    return(NA)
  }
  num<-as.numeric(as.character(substr(hhnum, 1, 2)))
  return(num)
}

mydf2$X2HHNUMBER<-as.numeric(lapply(mydf2$X2HHNUMBER, hh_fun))

mydf2$X2MTHID<-ifelse(mydf2$X2MTHID==-9, NA, mydf2$X2MTHID)
mydf2$X2MTHUTI<-ifelse(mydf2$X2MTHUTI==-9, NA, mydf2$X2MTHUTI)
mydf2$X2MTHEFF<-ifelse(mydf2$X2MTHEFF==-9, NA, mydf2$X2MTHEFF)
mydf2$X2MTHINT_R<-ifelse(mydf2$X2MTHINT_R==-9, NA, mydf2$X2MTHINT_R)

mydf2$X2STU30OCC_STEM1<-ifelse(mydf2$X2STU30OCC_STEM1=="Not a STEM occupation", 0, 
                               ifelse(mydf2$X2STU30OCC_STEM1=="Missing", NA, 1))

from<-c("Less than high school completion",                                              
        "Complete HS diploma/GED/alternative HS credential",                             
        "Start, but not complete certificate/diploma from school providing occ training",
        "Complete certificate/diploma from school providing occupational training",      
        "Start, but not complete Associate's degree",                                    
        "Complete Associate's degree",                                                   
        "Start, but not complete Bachelor's degree",                                     
        "Complete Bachelor's degree",                                                    
        "Start, but not complete Master's degree",                                       
        "Complete Master's degree",                                                      
        "Start, but not complete Ph.D./M.D./law degree/high level professional degree",  
        "Complete Ph.D./M.D./law degree/other high level professional degree",           
        "Don't know")

to<-c(1, 2, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 0)

mydf2<-mutate_at(mydf2, c("X2STUEDEXPCT"), mymap)

mydf2$no_expectation<-ifelse(mydf2$X2STUEDEXPCT==0, 1, 0)

mydf2$X2LOCALE<-ifelse(mydf2$X2LOCALE=="Missing", NA, mydf2$X2LOCALE)

mydf2$no_principal2<-ifelse(mydf2$X2SCHOOLCLI==-8, 1, 0)

mydf2$X2SCHOOLCLI<-ifelse(mydf2$X2SCHOOLCLI==-9, NA, mydf2$X2SCHOOLCLI)

mydf2$S2MCLUB<-ifelse(mydf2$S2MCLUB=="Yes", 1, 
                      ifelse(mydf2$S2MCLUB=="Missing", NA, 0))

mydf2$S2MCOMPETE<-ifelse(mydf2$S2MCOMPETE=="Yes", 1, 
                         ifelse(mydf2$S2MCOMPETE=="Missing", NA, 0))


from<-c("No time", "Less than 1/2 hour", "1/2 to 1 hour", "1 to 2 hours",            
        "2 to 3 hours", "4 to 6 hours", "7 to 9 hours", "More than 9 hours",       
        "Item legitimate skip/NA", "Missing")

to<-c(0, 0.5, 1, 2, 3, 6, 9, 10, 0, NA)

mydf2<-mutate_at(mydf2, c("S2MHOMEWRK"), mymap)


mydf2$S2MTUTORED<-ifelse(mydf2$S2MTUTORED=="Yes", 1, 
                         ifelse(mydf2$S2MTUTORED=="Missing", NA, 0))

mydf2$X2EVERDROP<-ifelse(mydf2$X2EVERDROP=="Yes", 1, 
                         ifelse(mydf2$X2EVERDROP=="Missing", NA, 0))

mydf2$X2BEHAVEIN<-ifelse(mydf2$X2BEHAVEIN==-9, NA, mydf2$X2BEHAVEIN)

mydf2$no_math2<-ifelse(mydf2$X2MEFFORT==-7, 1, 0)

mydf2$X2MEFFORT<-ifelse(mydf2$X2MEFFORT==-9, NA, mydf2$X2MEFFORT)
mydf2$X2PROBLEM<-ifelse(mydf2$X2PROBLEM==-9, NA, mydf2$X2PROBLEM)

from<-c("No",                                          
        "Yes",                                         
        "Item not administered: abbreviated interview",
        "Component not applicable",                    
        "Item legitimate skip/NA",                     
        "Unit non-response",                           
        "Missing")

to<-c(0, 1, 0, 0, 0, 0, NA)

mydf2<-mutate_at(mydf2, c("S2APMATH"), mymap)
mydf2<-mutate_at(mydf2, c("S2APSCIENCE"), mymap)

from<-c("No",                                          
        "Yes",                                         
        "Item not administered: abbreviated interview",
        "Component not applicable",                    
        "Item legitimate skip/NA",                     
        "Unit non-response",                           
        "Missing" )

to<-c(0, 1, 0, 0, 0, 0, NA)

mydf2<-mutate_at(mydf2, c("S2MSPR12"), mymap)


mydf2$S2HSJOBEVER<-ifelse(mydf2$S2HSJOBEVER=="Missing", NA, 
                          ifelse(mydf2$S2HSJOBEVER=="Yes", 1, 0))

mydf2$S2HSJOBHR<-ifelse(mydf2$S2HSJOBHR==-9, NA, ifelse(mydf2$S2HSJOBHR==-7, 0, mydf2$S2HSJOBHR))



from<-c("Pre-Algebra",                                                       
        "Algebra I, 1A or 1B",                                               
        "Algebra II",                                                        
        "Algebra III",                                                       
        "Geometry",                                                          
        "Analytic Geometry",                                                 
        "Trigonometry",                                                      
        "Pre-calculus or Analysis and Functions",                            
        "AP Calculus AB or BC",                                              
        "Other Calculus",                                                    
        "AP Statistics or Probability",                                      
        "Other Statistics or Probability",                                   
        "Integrated Math I",                                                 
        "Integrated Math II",                                                
        "Integrated Math III or above",                                      
        "IB mathematics standard or higher level",                           
        "Business/Consumer/General/Applied/Technical/Functional/Review math",
        "Other math course",                                                 
        "Component not applicable",                                          
        "Item legitimate skip/NA",                                           
        "Unit non-response",                                                 
        "Missing")

to<-c(1,2,3,4,5,6,7,8,9,8,8,6,2,3,4,9,3,4,0,0,0,NA)

mydf2<-mutate_at(mydf2, c("S2HIMATH12"), mymap)

mydf2$no_highest<-ifelse(is.na(mydf2$S2HIMATH12), 1, 0)

mydf2$moved_school<-ifelse(mydf2$X2NUMHS!="1 High school", 1, 0)




############################
############################
#Done Waves 1 and 2
############################
############################

mydf1$from<-1
mydf2$from<-2

my_combined<-merge(mydf1, mydf2, by="STU_ID", all.x=T)

my_combined_w1vals<-my_combined

my_combined_w2vals<-my_combined[!is.na(my_combined$X2SEX),]

X_W1<-rbind(my_combined_w1vals, my_combined_w2vals)
X_W2<-rbind(my_combined_w1vals, my_combined_w2vals)

X_W1<-X_W1[,c("X1SEX",
              "X1RACE",
              "X1DUALLANG",
              "X1STDOB",
              "X1TXMTH1.x",
              "X1PARPATTERN",
              "X1PAREDU",
              "X1HHNUMBER",
              "X1SES",
              "X1MTHID",
              "X1MTHUTI",
              "X1MTHEFF",
              "X1MTHINT",
              "X1SCHOOLBEL",
              "X1SCHOOLENG",
              "X1STU30OCC_STEM1",
              "X1STUEDEXPCT",
              "X1CONTROL",
              "X1LOCALE",
              "X1REGION",
              "X1SCHOOLCLI",
              "S1MCLUB",
              "S1MCOMPETE",
              "S1MCAMP",
              "S1MTUTOR",
              "S1MFALL09",
              "S1HRMHOMEWK",
              "S1HRWORK",
              "X2TXMTH1.y",
              "STU_ID",
              "S1M8",
              "S1M8GRADE",
              "S1S8GRADE",
              "parent_missing",
              "no_principal")]


X_W2<-X_W2[,c("X2SEX",
              "X2RACE",           
              "X2DUALLANG",       
              "X2STDOB",          
              "X2TXMTH1.y",         
              "X2PARPATTERN",     
              "X2PAREDU",         
              "X2HHNUMBER",       
              "X2SES",            
              "X2MTHID",         
              "X2MTHINT_R",       
              "X2MTHUTI",         
              "X2MTHEFF",         
              "X2STU30OCC_STEM1", 
              "X2STUEDEXPCT",    
              "X2CONTROL",        
              "X2LOCALE",         
              "X2REGION",         
              "X2SCHOOLCLI",      
              "S2MCLUB",         
              "S2MCOMPETE",       
              "S2MHOMEWRK",       
              "S2HSJOBHR",        
              "S2MTUTORED",       
              "X2EVERDROP",      
              "X2BEHAVEIN",       
              "X2MEFFORT",        
              "X2PROBLEM",        
              "S2APMATH",         
              "S2APSCIENCE",     
              "S2MSPR12",         
              "S2HSJOBEVER",      
              "X1TXMTH1.x",          
              "STU_ID",           
              "S2HIMATH12",      
              "X2NUMHS",          
              "no_expectation",   
              "no_principal2",    
              "no_math2",         
              "no_highest",      
              "moved_school",
              
              "X1HHNUMBER",
              "X1SES",
              "X1MTHID",
              "X1MTHUTI",
              "X1MTHEFF",
              "X1MTHINT",
              "X1SCHOOLBEL",
              "X1SCHOOLENG",
              "X1STU30OCC_STEM1",
              "X1STUEDEXPCT",
              "X1CONTROL",
              "X1LOCALE",
              "X1REGION",
              "X1SCHOOLCLI",
              "S1MCLUB",
              "S1MCOMPETE",
              "S1MCAMP",
              "S1MTUTOR",
              "S1MFALL09",
              "S1HRMHOMEWK",
              "S1HRWORK",
              "S1M8",
              "S1M8GRADE",
              "S1S8GRADE",
              "parent_missing",
              "no_principal")]

hist(X_W1$X1TXMTH1.x[1:21444])
hist(X_W2$X2TXMTH1.y[21445:40067])



#########################
#Modelling
#########################

y0<-my_combined_w1vals$X1TXMTH1.x
y1<-my_combined_w2vals$X2TXMTH1.y
y_vals<-c(y0, y1)


n_iter<-5000
n_burn<-3000
keep_every<-2
n_tree_mu<-100
n_tree_growth<-70
n_tree_treat<-30

Z_growth<-c(rep(0, 21444), rep(1, 18623))
Z_treat<-c(rep(0, 21444), my_combined_w2vals$S2HSJOBHR>=20)

X_W1<-X_W1[,-c(5,
               18,
               19,
               20,
               21,
               29,
               30,
               35)]

X_W2<-X_W2[,-c(5,
               16,
               17,
               18,
               19,
               23,
               28,
               32,
               34,
               38,
               
               29,
               30,
               31,
               35,
               40,
               52,
               53,
               54,
               55,
               67)]

X_W1<-as.matrix(one_hot(as.data.table(X_W1)))
X_W2<-as.matrix(one_hot(as.data.table(X_W2)))

valid_z_treat<-!is.na(Z_treat[21445:40067])

p_mod<-bartMachine(data.frame(X_W2[c(21445:40067)[valid_z_treat],]), as.factor(Z_treat[c(21445:40067)[valid_z_treat]]), num_trees=50, run_in_sample = F, mem_cache_for_speed = F, num_burn_in = 500, num_iterations_after_burn_in = 500, use_missing_data = T)
p<-predict(p_mod, data.frame(X_W2[c(21445:40067),]), type="prob")
p_scores<-rep(mean(p), 40067)
p<-c(rep(0, 21444), p)

X_W2<-cbind(X_W2, p)


sourceCpp(file = "~/your_directory/LBCF_HSLS_Github.cpp")

z_missing<-is.na(Z_treat)

t<-Sys.time()
my_mod <- fast_lbcf2(X_W1, y_vals, Z_growth, Z_treat, X_W2, X_W2, 
                     0.95, 2,
                     0.95, 2,
                     0.25, 3,
                     n_tree_mu, n_tree_growth, n_tree_treat/(0.5^2),
                     3, 0.1, n_iter, 
                     n_tree_mu, n_tree_growth, n_tree_treat, 1,
                     p_scores,
                     z_missing,
                     n_burn,
                     keep_every)
print(Sys.time()-t)

save(my_mod, file = "LBCF_Paper_1_Results.RData", compress = "xz")

