require(tidyverse)
require(lubridate)
require(zoo)


### Reestruturação dos dados, os pacientes 773, 489 e 212 
#aparentam apresentar erros no registro da data de evolução
data <- readRDS('dados_att.rds') %>%group_by(pac,exame) %>% 
mutate(DT_EVOL=replace(DT_EVOL, pac==773, '2017-05-26')) %>% 
  mutate(DT_EVOL=replace(DT_EVOL, pac==489, '2017-04-14')) %>% 
  mutate(DT_EVOL=replace(DT_EVOL, pac==212, '2017-02-16'))%>% 
  mutate(DATA_FINAL =difftime(DT_EVOL,DT_IS,units='days') ) %>% 
  mutate(dt_sup=DIAS,dt_inic =lag(DIAS) ) %>% select(-DIAS) %>% 
  mutate(dt_ate_inter =difftime(DT_NOT,DT_IS,units='days') ) 


###### Estruturação dos dados
dados <- data %>% 
  select(pac,EVOLUCAO,exame,valor,SEXO,DT_NASC,FEBRE:DIARREIA,dt_inic,dt_sup,DATA_FINAL,OBITO,dt_ate_inter) %>% 
  mutate(Idade = year(as.Date('2019-12-30'))-year(DT_NASC)) %>% select(-DT_NASC) %>% 
  spread(exame,valor) %>% mutate(TARGET = ifelse(OBITO =='OBITO',1,0)) %>% select(-OBITO,-EVOLUCAO)

####### Seleção de variáveis
data_inic <- dados %>% group_by(pac) %>% filter(is.na(dt_inic)) %>% ungroup(pac) %>% 
  select(Idade,SEXO,TGO,TGP,HEMORRAGIA,DOR_ABDOM,MIALGIA,Leucocitos,Bastoes,Creatinina,Linfocitos,Ureia,DATA_FINAL,TARGET)

#require(writexl)
#write_xlsx(data_inic,'dados.xlsx')


dados2 <- dados %>% select(pac,dt_inic,dt_sup,DATA_FINAL,TGO,TGP,TARGET,Idade,Ht,Leucocitos,Bastoes,Creatinina,SEXO) %>% 
  mutate(DATA_FINAL = ifelse(DATA_FINAL<0,max(dt_sup),DATA_FINAL))  


dados_completos <- dados2 %>% group_by(pac)%>% filter(any(DATA_FINAL == dt_sup))

dados_faltantes <- dados2 %>% filter(!(pac %in% (dados_completos$pac %>% unique))) 

dados_faltantes_2 <- dados_faltantes %>% 
  group_by(pac) %>% 
  summarise(dt_inic = max(dt_sup),dt_sup =max(DATA_FINAL)) %>% 
  bind_rows(dados_faltantes,.) %>% arrange(pac) 

dados_final <- bind_rows(dados_completos,dados_faltantes_2)  %>% 
  group_by(pac) %>% mutate(TARGET = ifelse(dt_sup==max(dt_sup),max(TARGET,na.rm = TRUE),0))

dados_input <- dados_final %>% group_by(pac) %>% 
  mutate_at(.vars=vars(DATA_FINAL,TGO,TGP,Ht,Leucocitos,Idade,SEXO),.funs=list(~na.locf(., na.rm = FALSE))) %>% 
  mutate_at(.vars=vars(DATA_FINAL,TGO,TGP,Ht,Leucocitos,Idade,SEXO),.funs=list(~na.locf(., na.rm = FALSE,fromLast=TRUE)))

dados_input <- dados_input%>% 
  group_by(pac) %>% mutate(TARGET = ifelse(dt_sup==max(dt_sup),max(TARGET,na.rm = TRUE),0))

######### Frequência de óbitos
dados_input %>% select(pac,TARGET) %>% group_by(pac) %>% filter(TARGET==max(TARGET)) %>% unique %>% ungroup() %>% 
  select(TARGET) %>% table

############### Regressão auxiliar para discretizar as variáveis contínuas 
require(rpart)
require(rattle)
corte_TGO = rpart(TARGET~TGO, data=dados2,maxdepth=1)
fancyRpartPlot(corte_TGO)


corte_TGP = rpart(TARGET~TGP, data=dados2,maxdepth=1 )
fancyRpartPlot(corte_TGP)


corte_Idade= rpart(TARGET~Idade, data=dados2,maxdepth=1)
fancyRpartPlot(corte_Idade)


corte_Leucocitos= rpart(TARGET~Leucocitos, data=dados2,maxdepth=1 )
fancyRpartPlot(corte_Leucocitos)


corte_Ht= rpart(TARGET~Ht, data=dados2,maxdepth=1 )
fancyRpartPlot(corte_Ht)


corte_Creatinina= rpart(TARGET~Creatinina, data=dados2,maxdepth=1 )
fancyRpartPlot(corte_Creatinina)

########### Discretização das variáveis de acordo com o corte proposto pela
########## Árvore de decisão e truncamento do tempo de sobrevivência em 50 dias
dados_disc <- dados_input %>% filter(dt_sup!=0) %>% 
  mutate(Idade = ifelse(Idade<= 40,'<=40','>40')) %>% 
  mutate(dt_sup = ifelse(dt_sup <=50,dt_sup,50)) %>% 
  mutate(dt_inic = ifelse(is.na(dt_inic),0,dt_inic))
  
#  mutate(TGO = ifelse(TGO<= 3500,'<=3500','>3500')) %>% 
#  mutate(TGP = ifelse(TGP<= 3500,'<=3500','>3500')) %>% 
# mutate(Ht = ifelse(Ht >= 32.25,'>=32.25','<32.25')) 


require(survival)

######### Tempo de sobrevivência
surv_obj <- Surv(dados_disc$dt_inic,dados_disc$dt_sup, event=dados_disc$TARGET)

########## Ajuste do modelo de regressão de cox considerando algumas variáveis
fit <- coxph( surv_obj ~ TGO+TGP+Ht+Leucocitos+Creatinina+cluster(pac),
              data=dados_disc)
summary(fit)
ajuste_cox <- survfit(fit)
ajuste_cox$surv
ajuste_cox$time
############### Teste de hipótese da presença de proporcionalidade
############### H0 é a hipótese onde os dados são proporcionais, neste caso não se rejeita a hipótese
############### de proporcionalidade para nenhuma covariável
teste_proporcionalidade <- cox.zph(fit)
teste_proporcionalidade

############### Resíduos e Teste de schoenfeld para analisar a dependência dos coeficientes ao tempo,
############## Os testes indicam que não existem evidências para rejeitar-mos a hipótese de independência ao tempo.
require(survminer)
ggcoxzph(teste_proporcionalidade)


############ Resíduo de Martingale, deve ser observado para todas as variáveins individualmente
ggcoxfunctional(surv_obj ~ TGO,data = dados_disc)

ggcoxfunctional(surv_obj ~ TGP,data = dados_disc)



