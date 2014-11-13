window_familiar_response_15 <- read.csv('window-familiar-response-15.csv')
window_unfamiliar_response_15 <- read.csv('window-unfamiliar-response-15.csv')
window_familiar_response_19 <- read.csv('window-familiar-response-19.csv')
window_unfamiliar_response_19 <- read.csv('window-unfamiliar-response-19.csv')

# assign age groups
window_familiar_response_15$AgeGroup <- '15-months'
window_unfamiliar_response_15$AgeGroup <- '15-months'
window_familiar_response_19$AgeGroup <- '19-months'
window_unfamiliar_response_19$AgeGroup <- '19-months'

familiar <- rbind(window_familiar_response_15, window_familiar_response_19)
familiar$AgeGroup <- factor(familiar$AgeGroup)
familiar$Target_c <- -.5
familiar[which(familiar$Target == 'Animal'), 'Target_c'] <- .5
familiar$AgeGroup_c <- -.5
familiar[which(familiar$AgeGroup == '19-months'), 'AgeGroup_c'] <- .5

unfamiliar <- rbind(window_unfamiliar_response_15, window_unfamiliar_response_19)
unfamiliar$AgeGroup <- factor(unfamiliar$AgeGroup)
unfamiliar$Condition_c <- -.5
unfamiliar[which(unfamiliar$Condition == 'Informative'), 'Condition_c'] <- .5
unfamiliar$AgeGroup_c <- -.5
unfamiliar[which(unfamiliar$AgeGroup == '19-months'), 'AgeGroup_c'] <- .5

# familiar
library(lme4)
model <- lmer(ArcSin ~ Target_c*AgeGroup_c + (1 | Trial) + (1 | ParticipantName), data = familiar)
model_null <- lmer(ArcSin ~ Target_c + AgeGroup_c + (1 | Trial) + (1 | ParticipantName), data = familiar)
anova(model,model_null)

# unfamiliar
model <- lmer(ArcSin ~ Condition_c*AgeGroup_c + (1 | Trial) + (1 | ParticipantName), data = unfamiliar)
model_null <- lmer(ArcSin ~ Condition_c + AgeGroup_c + (1 | Trial) + (1 | ParticipantName), data = unfamiliar)
anova(model,model_null)