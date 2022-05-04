#----------------------------------------------
# life-cycle component - Survival probablity
#----------------------------------------------
    
# Population Numbers from Bell and Miller (2002)
pop = zeros(Max_Age) ; #  Max_Age in model is 81 so that starting at age 20, agents live until age 100
    pop(1)=	197316
    pop(2)=	197141
    pop(3)=	196959
    pop(4)=	196770
    pop(5)=	196580
    pop(6)=	196392
    pop(7)=	196205
    pop(8)=	196019
    pop(9)=	195830
    pop(10)=195634
    pop(11)=195429
    pop(12)=195211
    pop(13)=194982
    pop(14)=194739
    pop(15)=194482
    pop(16)=194211
    pop(17)=193924
    pop(18)=193619
    pop(19)=193294
    pop(20)=192945
    pop(21)=192571
    pop(22)=192169
    pop(23)=191736
    pop(24)=191271
    pop(25)=190774
    pop(26)=190243
    pop(27)=189673
    pop(28)=189060
    pop(29)=188402
    pop(30)=187699
    pop(31)=186944
    pop(32)=186133
    pop(33)=185258
    pop(34)=184313
    pop(35)=183290
    pop(36)=182181
    pop(37)=180976
    pop(38)=179665
    pop(39)=178238
    pop(40)=176689
    pop(41)=175009
    pop(42)=173187
    pop(43)=171214
    pop(44)=169064
    pop(45)=166714
    pop(46)=164147
    pop(47)=161343
    pop(48)=158304
    pop(49)=155048
    pop(50)=151604
    pop(51)=147990
    pop(52)=144189
    pop(53)=140180
    pop(54)=135960
    pop(55)=131532
    pop(56)=126888
    pop(57)=122012
    pop(58)=116888
    pop(59)=111506
    pop(60)=105861
    pop(61)=99957
    pop(62)=93806
    pop(63)=87434
    pop(64)=80882
    pop(65)=74204
    pop(66)=67462
    pop(67)=60721
    pop(68)=54053
    pop(69)=47533
    pop(70)=41241
    pop(71)=35259
    pop(72)=29663
    pop(73)=24522
    pop(74)=19890
    pop(75)=15805
    pop(76)=12284
    pop(77)=9331
    pop(78)=6924
    pop(79)=5016
    pop(80)=3550

    if (Max_Age>80)	
        pop(Max_Age)=2454
    end		

# Survival probabilities: surv(i)=prob(alive in i+1|alive in i)
    Surv_Pr = zeros(Max_Age) ;
    for (age=1:Max_Age-1) 
        Surv_Pr(age)= pop(age+1)/pop(age) ;
    end 
    Surv_Pr(Max_Age)=0 ;