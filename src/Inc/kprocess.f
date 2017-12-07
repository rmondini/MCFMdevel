      integer kcase
      common/kcase/kcase
      integer, parameter ::
     & kW_only=   1,
     & kW_1jet=   2,
     & kWbfrmc=   3,
     & kW_cjet=   4,
     & kWcjet0=   5,
     & kWbbmas=   6,
     & kWbbbar=   7,
     & kW_2jet=   8,
     & kW_3jet=   9,
     & kWbbjet=  10,
     & kZ_only=  11,
     & kttZbbl=  12,
     & kZ_1jet=  13,
     & kZ_2jet=  14,
     & kZ_3jet=  15,
     & kZbbmas=  16,
     & kZbbbar=  17,
     & kZbbjet=  18,
     & kWWqqbr=  19,
     & kWWnpol=  20
      integer, parameter ::
     & kWW_jet=  21,
     & kWZbbar=  22,
     & kZZlept=  23,
     & kZZ_jet=  24,
     & kWHbbar=  25,
     & kWHgaga=  26,
     & kWH__WW=  27,
     & kWH__ZZ=  28,
     & kZHbbar=  29,
     & kZHgaga=  30,
     & kZH__WW=  31,
     & kZH__ZZ=  32,
     & kggfus0=  33,
     & kHWW_4l=  34, 
     & kHWW_tb=  35, 
     & kHWWint=  36,
     & kHWWHpi=  37,
     & kggWW4l=  38,
     & kggWWbx=  39,
     & kHWW2lq=  40
      integer, parameter ::
     & kHWWdkW=  41,
     & kHZZ_4l=  42,
     & kHZZ_tb=  43,
     & kHZZint=  44,
     & kHZZHpi=  45,
     & kggZZ4l=  46,
     & kggZZbx=  47,
     & kHZZqgI=  48,
     & kHi_Zga=  49,
     & kHVV_tb=  50, 
     & kHVVint=  51,
     & kHVVHpi=  52,
     & kggVV4l=  53,
     & kggVVbx=  54,
     & kH_1jet=  55,
     & ktt_bbl=  56,
     & ktt_ldk=  57,
     & ktt_bbu=  58,
     & ktt_udk=  59,
     & ktt_bbh=  60
      integer, parameter ::
     & ktt_hdk=  61,
     & ktthWdk=  62,
     & kqq_ttg=  63,
     & ktt_tot=  64,
     & kbb_tot=  65,
     & kcc_tot=  66,
     & ktt_glu=  67,
     & kbq_tpq=  68,
     & kttdkay=  69,
     & kt_bbar=  70,
     & ktdecay=  71,
     & kW_tndk=  72,
     & kW_twdk=  73,
     & kWtdkay=  74,
     & kWtbwdk=  75,
     & kWtbndk=  76,
     & khttjet=  77,
     & kggfus1=  78,
     & kattjet=  79,
     & kHWWjet=  80
      integer, parameter ::
     & kHZZjet=  81,
     & kqq_Hqq=  82,
     & kqq_HWW=  83,
     & kqq_HZZ=  84,
     & kqq_Hgg=  85,
     & kqqHqqg=  86,
     & kqqZZqq=  87,
     & ktautau=  88,
     & kqqWWqq=  89,
     & kqqVVqq=  90,
     & kqqWWss=  91,
     & kqqWZqq=  92,
     & kqg_tbq=  93,
     & kqgtbqq=  94,
     & k4ftwdk=  95,
     & kdk_4ft=  96,
     & k4ftjet=  97,
     & kqq_tbg=  98,
     & kqqtbgg=  99,
     & kepem3j= 100
      integer, parameter ::
     & kWpWp2j= 101,
     & kWpWp3j= 102,
     & kWpmZjj= 103,
     & kWpmZbj= 104,
     & kWpmZbb= 105,
     & kgQ__ZQ= 106,
     & kgagajj= 107,
     & kggfus2= 108,
     & kHWW2jt= 109,
     & kHZZ2jt= 110,
     & kggfus3= 111,
     & kHWW3jt= 112,
     & kHZZ3jt= 113,
     & kdirgam= 114,
     & kgamjet= 115,
     & khflgam= 116,
     & kgamgam= 117,
     & kgmgmjt= 118,
     & ktrigam= 119,
     & kgmgmjj= 120
      integer, parameter ::
     & kfourga= 121,
     & ktwojet= 122,
     & kthrjet= 123,
     & kWgamma= 124,
     & kWgajet= 125,
     & kZgamma= 126,
     & kZgajet= 127,
     & kZ_2gam= 128,
     & kZ2gajt= 129,
     & kZga2jt= 130,
     & kW_bjet= 131,
     & kWcjetg= 132,
     & kZ_bjet= 133,
     & kZbjetg= 134,
     & kWcsbar= 135,
     & kWcs_ms= 136,
     & kW_2gam= 137,
     & kWbbjem= 138,
     & kWttmas= 139,
     & kqq_ttw= 140
      integer, parameter ::
     & kttwldk= 141,
     & kqq_ttz= 142,
     & kqqtthz= 143,
     & kH_tjet= 144,
     & kH_tdkj= 145,
     & kZ_tjet= 146,
     & kZt2jet= 147,
     & kZ_tdkj= 148,
     & kZtdk2j= 149,
     & kHHpair= 150,
     & kWH1jet= 151,
     & kZH1jet= 152,
     & ktottth= 153,
     & kqq_tth= 154,
     & ktth_ww= 155,
     & kdm_jet= 156,
     & kdm_gam= 157,
     & kdm2jet= 158,
     & kdm_gaj= 159,
     & kvlchk2= 160
      integer, parameter ::
     & kvlchk3= 161,
     & kvlchk4= 162,
     & kvlchk5= 163,
     & kvlchk6= 164,
     & kvlchk8= 165,
     & kvlchkm= 166,
     & kvlchm3= 167,
     & kvlchwt= 168,
     & kvlchwn= 169,
     & kvlchwg= 170,
     & kvlchwh= 171,
     & kWWqqdk= 172,
     & kHigaga= 173,
     & kdmgamj= 174,
     & kHgagaj= 175,
     & kHZZpjt= 176,
     & kHmZZ4l= 177,
     & kHZZ_jj= 178,
     & kW_cwdk= 179,
     & kZccmas= 180
      integer, parameter ::
     & ktotttz= 181,
     & kW_frag= 182,
     & kWb2jet= 183,
     & kWW2jet= 184,
     & kZ_frag= 185,
     & kqqttbb= 186,
     & kqqttgg= 187,
     & kqg_tbb= 188,
     & kHt2jet= 189,
     & kttbdkh= 190,
     & kttbdkl= 191,
     & khlljet= 192,
     & kZHlept= 193,
     & kWHbbdk= 194,
     & kZHbbdk= 195,
     & kgg2gam= 196,
     & kgam_2j= 197, 
     & kgam_3j= 198,
     & kZ2jetF= 199,
     & kHgagaI= 200,
     & kZHbbnn= 201,
     & kZHbbjt= 202,
     & kWHbbnn= 203,
     & kWHbbjt= 204,
     & kZHn2bc= 205,
     & kZHn2ac= 206
