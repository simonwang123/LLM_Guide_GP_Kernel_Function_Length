# 序贯大坝代理模型的物理特征进化

你在为三个独立的坝体变形约束响应提出衍生物理特征。当前模型将物理特征的 Matérn-5/2 ARD 核与原始八维设计的 Matérn-5/2 ARD 残差核相加。物理候选用于几何表达，模型的拟合和采用由本地程序处理。

每项响应的 `matched_pairs` 包含五个绝对留一预测残差最大的断面，以及从低残差样本中为各断面匹配的对照。带符号残差为有限元真值减去留一预测，正值表示低估。响应和残差以米计。几何原语、原始变量及其别名由上下文定义。

1. 分别对比每项响应的高残差组与低残差对照，寻找可重复出现的几何区别。根据残差符号描述方向。符号混合时，可以在力学解释中说明适用的断面子组。
2. 每项响应给出恰好两条机制不同的新公式。每条解释应指出支持它的观察与可能的力学机制，注明这些解释是候选假设。
3. 公式只可使用上下文中的原始变量、别名、几何原语、数值常量和 `+`、`-`、`*`、`/`、`**`、`sqrt`、`abs`、`log`、`log1p`、`min`、`max`。幂指数为绝对值不大于 3 的常量。遵守量纲一致与函数定义域，避免分母在几何范围内触零。公式应组合多个变化量或对几何原语无量纲化。已有特征值用于诊断，表达式应直接写成原始变量与原语的组合。
4. 返回一个纯 JSON 对象，键为下面 `responses` 中的三个响应全名。每个键对应两元素数组，每个元素只有 `name`、`expression`、`physical_rationale` 三个字符串字段。`name` 使用英文 snake_case，且与该响应的已有特征及另一候选不同。公式只依赖设计几何，真值、预测、残差和样本编号用于诊断。

## 当前已观测前缀

```json
{
  "sample_count": 32,
  "response_unit": "m",
  "residual_definition": "observed minus LOO prediction",
  "original_variable_units": {
    "upstream_slope": "1",
    "downstream_slope": "1",
    "crest_width": "m",
    "modulus_zone_height": "m",
    "secondary_top_width": "m",
    "secondary_height": "m",
    "secondary_bottom_width": "m",
    "secondary_right_slope_distance": "m"
  },
  "aliases": {
    "x1": "upstream_slope",
    "x2": "downstream_slope",
    "x3": "crest_width",
    "x4": "modulus_zone_height",
    "x5": "secondary_top_width",
    "x6": "secondary_height",
    "x7": "secondary_bottom_width",
    "x8": "secondary_right_slope_distance",
    "d_ds_perp": "secondary_right_slope_distance"
  },
  "primitive_catalog": {
    "A_tot": {
      "unit": "m^2",
      "desc": "坝体主体域面积（GeoStudio 材料 4、5、6、12 的并集）"
    },
    "A_mod": {
      "unit": "m^2",
      "desc": "坝顶增模区(Material 12)面积（顶部高模量加强帽盖）"
    },
    "A_sec": {
      "unit": "m^2",
      "desc": "次堆石区(Material 5)面积（内部低模量软弱压缩区）"
    },
    "A_drain": {
      "unit": "m^2",
      "desc": "底部排水区(Material 6)面积（底部高刚度支承区）"
    },
    "A_main": {
      "unit": "m^2",
      "desc": "主堆石区(Material 4)面积（高刚度主承重与抗推骨架）"
    },
    "A_ds_shell": {
      "unit": "m^2",
      "desc": "次堆石区顶底高程之间、其右腰下游侧的材料 4 面积"
    },
    "B_bot": {
      "unit": "m",
      "desc": "坝体主体域在坝底高程的水平宽度"
    },
    "B_mod_bot": {
      "unit": "m",
      "desc": "坝顶增模区底面宽度"
    },
    "h_drain": {
      "unit": "m",
      "desc": "底部排水区高度（次堆石区底面至坝基的距离）"
    },
    "d_ds_perp": {
      "unit": "m",
      "desc": "次堆石区右腰到 GeoStudio 下游参考线的法向距离，等于 x8"
    },
    "d_us_min": {
      "unit": "m",
      "desc": "次堆石区左腰距主体域上游边界的最小水平净宽"
    },
    "X_c_sec": {
      "unit": "m",
      "desc": "次堆石区形心相对坝轴线的横坐标，正值偏下游"
    },
    "Y_c_sec": {
      "unit": "m",
      "desc": "次堆石区形心相对坝底的高程"
    },
    "W_sec_mean": {
      "unit": "m",
      "desc": "次堆石区平均水平宽度 0.5*(x5+x7)"
    },
    "H_dam": {
      "unit": "m",
      "desc": "大坝总高度常量"
    }
  },
  "responses": {
    "crest_max_settlement_m": {
      "current_formulas": [
        {
          "name": "narrow_crest_broad_base_index",
          "expression": "B_bot / crest_width",
          "physical_rationale": "候选假设：高残差断面相对低残差对照多数表现为坝顶更窄、坝底更宽，即更大的坝底宽度与坝顶宽度之比。正残差子组中这种几何可能使坝顶竖向荷载向更宽、更深的坝体扩散，放大坝顶下方压缩区而让模型低估沉降；负残差子组也出现该几何，说明方向会被帽盖和壳体刚度调节，单独使用时主要表征荷载扩散尺度。"
        },
        {
          "name": "soft_area_ratio",
          "expression": "A_sec / A_tot",
          "physical_rationale": "低模量次堆石占比表征主体断面的整体可压缩性。"
        },
        {
          "name": "stiff_support_to_soft_area",
          "expression": "(A_mod + A_ds_shell) / A_sec",
          "physical_rationale": "候选假设：高残差断面相对对照通常降低 A_mod 与 A_ds_shell 之和相对于 A_sec 的比例，正残差子组尤其明显；负残差中帽盖显著增大的断面可能反向。力学上，顶部增模区和下游壳体构成绕过软区的刚度骨架，骨架面积不足会使软区承担更大压缩应变，导致观测沉降高于预测。"
        }
      ],
      "loo_nrmse": 0.3128976760475957,
      "loo_formula_fallback": null,
      "top5_signs": {
        "underestimated": 3,
        "overestimated": 2,
        "zero": 0
      },
      "top5_direction": "mixed_or_tied",
      "matched_pairs": [
        {
          "high_residual": {
            "sample_index": 7,
            "design_id": "lhs_007",
            "original_variables": {
              "upstream_slope": 1.3964213037264912,
              "downstream_slope": 1.55273578190388,
              "crest_width": 9.43626151664104,
              "modulus_zone_height": 12.050513005797814,
              "secondary_top_width": 12.720244476030498,
              "secondary_height": 20.88160543608512,
              "secondary_bottom_width": 35.83143862047909,
              "secondary_right_slope_distance": 13.627449547578786
            },
            "geometric_primitives": {
              "A_tot": 8042.333623728347,
              "A_mod": 294.3714524013453,
              "A_sec": 506.91854484567375,
              "A_drain": 3132.356073842461,
              "A_main": 4108.687552638867,
              "A_ds_shell": 522.9673958031354,
              "B_bot": 220.27388134400002,
              "B_mod_bot": 42.58935197700001,
              "h_drain": 38.06788155800001,
              "d_ds_perp": 13.627449547578786,
              "d_us_min": 4.945136433000016,
              "X_c_sec": 8.215963533000007,
              "Y_c_sec": 46.85203440499998,
              "W_sec_mean": 24.275841548499997,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "narrow_crest_broad_base_index": 23.34334216527833,
              "soft_area_ratio": 0.06303127531914937,
              "stiff_support_to_soft_area": 1.612367226480758
            },
            "loo_prediction_m": 0.08504854869133893,
            "observed_m": 0.07708631363929866,
            "observed_minus_loo_m": -0.007962235052040273
          },
          "matched_low_residual": {
            "sample_index": 32,
            "design_id": "eval_0032",
            "original_variables": {
              "upstream_slope": 1.3320104010403158,
              "downstream_slope": 1.3892031960189342,
              "crest_width": 7.83259486220777,
              "modulus_zone_height": 17.462718952447176,
              "secondary_top_width": 18.06828137487173,
              "secondary_height": 24.233660977333784,
              "secondary_bottom_width": 33.185473438352346,
              "secondary_right_slope_distance": 10.66262403037399
            },
            "geometric_primitives": {
              "A_tot": 7306.114505716314,
              "A_mod": 494.31824134870135,
              "A_sec": 621.0330589964988,
              "A_drain": 1779.5526228276067,
              "A_main": 4411.210582543508,
              "A_ds_shell": 438.5025185213277,
              "B_bot": 201.80239653600003,
              "B_mod_bot": 52.634710216,
              "h_drain": 29.303620070000022,
              "d_ds_perp": 10.66262403037399,
              "d_us_min": 16.62668582699999,
              "X_c_sec": 19.091129449999983,
              "Y_c_sec": 40.229172285000004,
              "W_sec_mean": 25.6268774065,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "narrow_crest_broad_base_index": 25.76443695686286,
              "soft_area_ratio": 0.08500182395315617,
              "stiff_support_to_soft_area": 1.5020468658743142
            },
            "loo_prediction_m": 0.07262432573946621,
            "observed_m": 0.07437356009175546,
            "observed_minus_loo_m": 0.0017492343522892417
          }
        },
        {
          "high_residual": {
            "sample_index": 25,
            "design_id": "eval_0025",
            "original_variables": {
              "upstream_slope": 1.6882901702076196,
              "downstream_slope": 1.6950646217912435,
              "crest_width": 11.977469097822905,
              "modulus_zone_height": 18.36188363377005,
              "secondary_top_width": 11.75962993875146,
              "secondary_height": 27.23224811255932,
              "secondary_bottom_width": 25.882059396244586,
              "secondary_right_slope_distance": 13.828325916081667
            },
            "geometric_primitives": {
              "A_tot": 9408.289879923997,
              "A_mod": 766.9773790584494,
              "A_sec": 512.5339116636007,
              "A_drain": 1520.5571967507958,
              "A_main": 6608.221392451151,
              "A_ds_shell": 737.2519028341619,
              "B_bot": 254.94571919599997,
              "B_mod_bot": 73.42945398900002,
              "h_drain": 25.405868253999984,
              "d_ds_perp": 13.828325916081667,
              "d_us_min": 34.727601421,
              "X_c_sec": 30.76631473499998,
              "Y_c_sec": 37.31915736899998,
              "W_sec_mean": 18.820844667500026,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "narrow_crest_broad_base_index": 21.285441616571603,
              "soft_area_ratio": 0.054476840978005786,
              "stiff_support_to_soft_area": 2.934887326791959
            },
            "loo_prediction_m": 0.05048949701133716,
            "observed_m": 0.05785997394491101,
            "observed_minus_loo_m": 0.007370476933573851
          },
          "matched_low_residual": {
            "sample_index": 23,
            "design_id": "eval_0023",
            "original_variables": {
              "upstream_slope": 1.3465167950838806,
              "downstream_slope": 1.3504548892378807,
              "crest_width": 9.308618433773518,
              "modulus_zone_height": 18.278893544338644,
              "secondary_top_width": 29.175809305161238,
              "secondary_height": 24.055996481329203,
              "secondary_bottom_width": 27.749729226343334,
              "secondary_right_slope_distance": 3.510688710026443
            },
            "geometric_primitives": {
              "A_tot": 7344.622640850989,
              "A_mod": 559.8062705784446,
              "A_sec": 684.7002773229959,
              "A_drain": 944.4711616301922,
              "A_main": 5155.644931319356,
              "A_ds_shell": 138.2280912578236,
              "B_bot": 201.48451856399998,
              "B_mod_bot": 55.86894898700001,
              "h_drain": 28.665109974000018,
              "d_ds_perp": 3.510688710026443,
              "d_us_min": 21.101340232000013,
              "X_c_sec": 29.066374193,
              "Y_c_sec": 40.793548714999986,
              "W_sec_mean": 28.4627692655,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "narrow_crest_broad_base_index": 21.644943339064593,
              "soft_area_ratio": 0.09322470476763156,
              "stiff_support_to_soft_area": 1.0194743378305688
            },
            "loo_prediction_m": 0.07119354987691902,
            "observed_m": 0.07307196200591454,
            "observed_minus_loo_m": 0.0018784121289955186
          }
        },
        {
          "high_residual": {
            "sample_index": 19,
            "design_id": "lhs_019",
            "original_variables": {
              "upstream_slope": 1.6745232365833786,
              "downstream_slope": 1.4008319126769921,
              "crest_width": 9.165674889626958,
              "modulus_zone_height": 16.760990775170942,
              "secondary_top_width": 24.924136868277312,
              "secondary_height": 36.21510900519744,
              "secondary_bottom_width": 57.11954572428753,
              "secondary_right_slope_distance": 10.172838528352186
            },
            "geometric_primitives": {
              "A_tot": 8367.75443834113,
              "A_mod": 547.9646538156717,
              "A_sec": 1485.6104541145198,
              "A_drain": 1491.4065154909613,
              "A_main": 4842.7728149199775,
              "A_ds_shell": 626.5786803350254,
              "B_bot": 229.341951425,
              "B_mod_bot": 59.01079680599997,
              "h_drain": 18.023900219999973,
              "d_ds_perp": 10.172838528352186,
              "d_us_min": 16.988913295000003,
              "X_c_sec": 21.663753582000027,
              "Y_c_sec": 33.76287860899998,
              "W_sec_mean": 41.02184129599998,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "narrow_crest_broad_base_index": 25.02182918189172,
              "soft_area_ratio": 0.17753992006594252,
              "stiff_support_to_soft_area": 0.790613266686232
            },
            "loo_prediction_m": 0.07707568716507636,
            "observed_m": 0.08383617259362591,
            "observed_minus_loo_m": 0.006760485428549554
          },
          "matched_low_residual": {
            "sample_index": 20,
            "design_id": "lhs_020",
            "original_variables": {
              "upstream_slope": 1.4420337151416434,
              "downstream_slope": 1.6952251044089008,
              "crest_width": 10.899606022680448,
              "modulus_zone_height": 13.154257709251336,
              "secondary_top_width": 23.49373414724033,
              "secondary_height": 32.683541574470254,
              "secondary_bottom_width": 52.227206782238454,
              "secondary_right_slope_distance": 7.27158587734832
            },
            "geometric_primitives": {
              "A_tot": 8659.800493040262,
              "A_mod": 386.1615117911777,
              "A_sec": 1237.4142604653896,
              "A_drain": 1956.4091925520615,
              "A_main": 5079.815528231634,
              "A_ds_shell": 461.85458236732734,
              "B_bot": 235.65675416000002,
              "B_mod_bot": 50.418226704000006,
              "h_drain": 25.162200715999973,
              "d_ds_perp": 7.27158587734832,
              "d_us_min": 12.972606106,
              "X_c_sec": 29.586810307000007,
              "Y_c_sec": 39.436926298,
              "W_sec_mean": 37.8604704645,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "narrow_crest_broad_base_index": 21.62066717545878,
              "soft_area_ratio": 0.14289177463844333,
              "stiff_support_to_soft_area": 0.6853130121836218
            },
            "loo_prediction_m": 0.07927347101180572,
            "observed_m": 0.07761399443976472,
            "observed_minus_loo_m": -0.0016594765720410048
          }
        },
        {
          "high_residual": {
            "sample_index": 13,
            "design_id": "lhs_013",
            "original_variables": {
              "upstream_slope": 1.6446159216739957,
              "downstream_slope": 1.6795317296966799,
              "crest_width": 7.589823503377536,
              "modulus_zone_height": 11.310364914746476,
              "secondary_top_width": 27.938447008999383,
              "secondary_height": 24.297527356364657,
              "secondary_bottom_width": 65.86501425442128,
              "secondary_right_slope_distance": 8.027170857060547
            },
            "geometric_primitives": {
              "A_tot": 8935.356992263625,
              "A_mod": 279.49980368069726,
              "A_sec": 1139.5960830739239,
              "A_drain": 3871.8024894206746,
              "A_main": 3644.4586160883296,
              "A_ds_shell": 378.5524835054062,
              "B_bot": 246.17674519499997,
              "B_mod_bot": 43.878225472000025,
              "h_drain": 35.392107728999974,
              "d_ds_perp": 8.027170857060547,
              "d_us_min": 0.49380531299999575,
              "X_c_sec": 10.12346508600001,
              "Y_c_sec": 45.90354411200002,
              "W_sec_mean": 46.90173063150003,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "narrow_crest_broad_base_index": 32.43510802134585,
              "soft_area_ratio": 0.12753783470102026,
              "stiff_support_to_soft_area": 0.5774434441816317
            },
            "loo_prediction_m": 0.08443318389080556,
            "observed_m": 0.07951806597548103,
            "observed_minus_loo_m": -0.004915117915324535
          },
          "matched_low_residual": {
            "sample_index": 29,
            "design_id": "eval_0029",
            "original_variables": {
              "upstream_slope": 1.3082879073917866,
              "downstream_slope": 1.548298779502511,
              "crest_width": 6.455660009756684,
              "modulus_zone_height": 15.85861544124782,
              "secondary_top_width": 20.622952468693256,
              "secondary_height": 37.372649908065796,
              "secondary_bottom_width": 27.190962894819677,
              "secondary_right_slope_distance": 4.046557015739381
            },
            "geometric_primitives": {
              "A_tot": 7578.107427309267,
              "A_mod": 415.0985016513749,
              "A_sec": 893.4663598054761,
              "A_drain": 645.1210968156714,
              "A_main": 5624.421469036744,
              "A_ds_shell": 270.7107643710405,
              "B_bot": 210.44307032999998,
              "B_mod_bot": 49.34131100000002,
              "h_drain": 17.76873465099999,
              "d_ds_perp": 4.046557015739381,
              "d_us_min": 21.701797064000004,
              "X_c_sec": 43.00476230600003,
              "Y_c_sec": 35.599437163000005,
              "W_sec_mean": 23.90695768149999,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "narrow_crest_broad_base_index": 32.59822698406505,
              "soft_area_ratio": 0.11790098891785125,
              "stiff_support_to_soft_area": 0.7675826386700542
            },
            "loo_prediction_m": 0.06722706965000042,
            "observed_m": 0.06854139129052197,
            "observed_minus_loo_m": 0.0013143216405215596
          }
        },
        {
          "high_residual": {
            "sample_index": 2,
            "design_id": "lhs_002",
            "original_variables": {
              "upstream_slope": 1.569363259687807,
              "downstream_slope": 1.4665851746995955,
              "crest_width": 6.669684274972619,
              "modulus_zone_height": 7.144992267623804,
              "secondary_top_width": 15.054269429416003,
              "secondary_height": 22.4629175393885,
              "secondary_bottom_width": 51.30574787263373,
              "secondary_right_slope_distance": 5.155855220570123
            },
            "geometric_primitives": {
              "A_tot": 8083.082706026641,
              "A_mod": 106.02391874393064,
              "A_sec": 745.3197982798283,
              "A_drain": 3894.8369460715667,
              "A_main": 3336.9020429313186,
              "A_ds_shell": 202.66168602689874,
              "B_bot": 223.92986391,
              "B_mod_bot": 25.91713976700001,
              "h_drain": 41.392090193,
              "d_ds_perp": 5.155855220570123,
              "d_us_min": 1.969850661999999,
              "X_c_sec": 9.817879002000012,
              "Y_c_sec": 50.578356681,
              "W_sec_mean": 33.18000865150002,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "narrow_crest_broad_base_index": 33.57428248144764,
              "soft_area_ratio": 0.09220737005748161,
              "stiff_support_to_soft_area": 0.4141653092850409
            },
            "loo_prediction_m": 0.08555929445946754,
            "observed_m": 0.09035779592611756,
            "observed_minus_loo_m": 0.004798501466650021
          },
          "matched_low_residual": {
            "sample_index": 24,
            "design_id": "eval_0024",
            "original_variables": {
              "upstream_slope": 1.303603621944785,
              "downstream_slope": 1.403586372360587,
              "crest_width": 6.122953802347183,
              "modulus_zone_height": 17.225281218998134,
              "secondary_top_width": 13.17300047725439,
              "secondary_height": 24.71127737313509,
              "secondary_bottom_width": 43.46261053811759,
              "secondary_right_slope_distance": 3.6680213138461113
            },
            "geometric_primitives": {
              "A_tot": 7146.540215651556,
              "A_mod": 449.64511115419907,
              "A_sec": 699.7691464892857,
              "A_drain": 1968.8558545715432,
              "A_main": 4028.2701034365264,
              "A_ds_shell": 152.21815763154154,
              "B_bot": 199.055008872,
              "B_mod_bot": 49.97987992299997,
              "h_drain": 29.063441408000017,
              "d_ds_perp": 3.6680213138461113,
              "d_us_min": 30.80764248699998,
              "X_c_sec": 29.75746385299999,
              "Y_c_sec": 39.216417192999984,
              "W_sec_mean": 28.317805507999992,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "narrow_crest_broad_base_index": 32.509637553641824,
              "soft_area_ratio": 0.09791719144834998,
              "stiff_support_to_soft_area": 0.8600883188481014
            },
            "loo_prediction_m": 0.06766510251600441,
            "observed_m": 0.06757443292644877,
            "observed_minus_loo_m": -9.066958955564197e-05
          }
        }
      ]
    },
    "dam_max_settlement_m": {
      "current_formulas": [
        {
          "name": "soft_core_eccentricity",
          "expression": "X_c_sec / B_bot",
          "physical_rationale": "软区相对坝轴线的有符号偏心反映骨架承载分布的不对称性。"
        },
        {
          "name": "deep_soft_drainage_deficit",
          "expression": "((H_dam - Y_c_sec) / H_dam) * (A_sec / (A_sec + A_drain))",
          "physical_rationale": "候选假设：正残差断面的次堆石形心更低且 A_sec/(A_sec+A_drain) 更大，而多数负残差对照形心更高或底部排水支承更强；个别负残差断面因主骨架和帽盖很大而不服从该方向。力学上，低位软区若缺少底部高刚度排水支承，会形成更长的压缩柱并增大整体沉降，使模型低估。"
        }
      ],
      "loo_nrmse": 0.2166886228361784,
      "loo_formula_fallback": null,
      "top5_signs": {
        "underestimated": 2,
        "overestimated": 3,
        "zero": 0
      },
      "top5_direction": "mixed_or_tied",
      "matched_pairs": [
        {
          "high_residual": {
            "sample_index": 10,
            "design_id": "lhs_010",
            "original_variables": {
              "upstream_slope": 1.5468086093385431,
              "downstream_slope": 1.3012624657189515,
              "crest_width": 11.291715983318776,
              "modulus_zone_height": 13.534298401099962,
              "secondary_top_width": 26.23483907948985,
              "secondary_height": 27.391521289178613,
              "secondary_bottom_width": 63.350919984777626,
              "secondary_right_slope_distance": 8.592072960145197
            },
            "geometric_primitives": {
              "A_tot": 7897.915705416024,
              "A_mod": 372.63251196399375,
              "A_sec": 1226.9451132903607,
              "A_drain": 2947.1016370375714,
              "A_main": 3351.236443124097,
              "A_ds_shell": 381.73199744837643,
              "B_bot": 214.648971028,
              "B_mod_bot": 47.245817955999996,
              "h_drain": 30.074180309999974,
              "d_ds_perp": 8.592072960145197,
              "d_us_min": 7.236450769000015,
              "X_c_sec": 9.331825449999997,
              "Y_c_sec": 41.878520691000006,
              "W_sec_mean": 44.792879532000015,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_core_eccentricity": 0.043474820332507914,
              "deep_soft_drainage_deficit": 0.12759440893394416
            },
            "loo_prediction_m": 0.18821979450643495,
            "observed_m": 0.17887868717138825,
            "observed_minus_loo_m": -0.0093411073350467
          },
          "matched_low_residual": {
            "sample_index": 1,
            "design_id": "lhs_001",
            "original_variables": {
              "upstream_slope": 1.6094092387020262,
              "downstream_slope": 1.398042551563892,
              "crest_width": 10.294596895348278,
              "modulus_zone_height": 17.22436933287436,
              "secondary_top_width": 22.75627311159097,
              "secondary_height": 28.75351770691302,
              "secondary_bottom_width": 46.6598200010875,
              "secondary_right_slope_distance": 12.621874746349105
            },
            "geometric_primitives": {
              "A_tot": 8262.411076914119,
              "A_mod": 581.5095874314916,
              "A_sec": 997.9784312431392,
              "A_drain": 1974.2684502131478,
              "A_main": 4708.654608026342,
              "A_ds_shell": 619.0967573010715,
              "B_bot": 225.44602486499997,
              "B_mod_bot": 60.22160053099998,
              "h_drain": 25.022112960000015,
              "d_ds_perp": 12.621874746349105,
              "d_us_min": 16.09770578199999,
              "X_c_sec": 15.802945496999996,
              "Y_c_sec": 37.74865155499998,
              "W_sec_mean": 34.70804655650002,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_core_eccentricity": 0.07009635901304095,
              "deep_soft_drainage_deficit": 0.16448591730430925
            },
            "loo_prediction_m": 0.16958088876560878,
            "observed_m": 0.1688936901258416,
            "observed_minus_loo_m": -0.0006871986397671836
          }
        },
        {
          "high_residual": {
            "sample_index": 11,
            "design_id": "lhs_011",
            "original_variables": {
              "upstream_slope": 1.3717017078040479,
              "downstream_slope": 1.320912292949585,
              "crest_width": 9.763338568555401,
              "modulus_zone_height": 5.227712698589564,
              "secondary_top_width": 14.210190616724496,
              "secondary_height": 29.327369356524045,
              "secondary_bottom_width": 49.60905740491252,
              "secondary_right_slope_distance": 2.713831948671564
            },
            "geometric_primitives": {
              "A_tot": 7364.980022583283,
              "A_mod": 68.12974081013377,
              "A_sec": 935.825329374874,
              "A_drain": 2780.5265482474133,
              "A_main": 3580.498404150864,
              "A_ds_shell": 126.00460934764087,
              "B_bot": 201.616770115,
              "B_mod_bot": 20.240492020000005,
              "h_drain": 36.444917944999986,
              "d_ds_perp": 2.713831948671564,
              "d_us_min": 1.906027913999992,
              "X_c_sec": 16.229618721999998,
              "Y_c_sec": 48.39741523800001,
              "W_sec_mean": 31.909624011000005,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_core_eccentricity": 0.08049736494014263,
              "deep_soft_drainage_deficit": 0.08712244618321494
            },
            "loo_prediction_m": 0.16282413905793405,
            "observed_m": 0.1719201429519541,
            "observed_minus_loo_m": 0.009096003894020044
          },
          "matched_low_residual": {
            "sample_index": 28,
            "design_id": "eval_0028",
            "original_variables": {
              "upstream_slope": 1.3775937989354134,
              "downstream_slope": 1.3864084053784609,
              "crest_width": 6.271950423717499,
              "modulus_zone_height": 13.629157315008342,
              "secondary_top_width": 13.016253504902124,
              "secondary_height": 26.466528959572315,
              "secondary_bottom_width": 31.796314222738147,
              "secondary_right_slope_distance": 11.256418767385185
            },
            "geometric_primitives": {
              "A_tot": 7312.249949311609,
              "A_mod": 297.46670565444634,
              "A_sec": 593.0165607714652,
              "A_drain": 1921.08893760141,
              "A_main": 4500.677745284286,
              "A_ds_shell": 504.6623068138584,
              "B_bot": 203.40810903400003,
              "B_mod_bot": 41.104390947999974,
              "h_drain": 30.90431372500001,
              "d_ds_perp": 11.256418767385185,
              "d_us_min": 9.184546398000009,
              "X_c_sec": 15.678812819000001,
              "Y_c_sec": 42.288978411000016,
              "W_sec_mean": 22.406283863500008,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_core_eccentricity": 0.07708056917425676,
              "deep_soft_drainage_deficit": 0.10107920965855705
            },
            "loo_prediction_m": 0.15921818754603995,
            "observed_m": 0.1602549082414876,
            "observed_minus_loo_m": 0.0010367206954476627
          }
        },
        {
          "high_residual": {
            "sample_index": 23,
            "design_id": "eval_0023",
            "original_variables": {
              "upstream_slope": 1.3465167950838806,
              "downstream_slope": 1.3504548892378807,
              "crest_width": 9.308618433773518,
              "modulus_zone_height": 18.278893544338644,
              "secondary_top_width": 29.175809305161238,
              "secondary_height": 24.055996481329203,
              "secondary_bottom_width": 27.749729226343334,
              "secondary_right_slope_distance": 3.510688710026443
            },
            "geometric_primitives": {
              "A_tot": 7344.622640850989,
              "A_mod": 559.8062705784446,
              "A_sec": 684.7002773229959,
              "A_drain": 944.4711616301922,
              "A_main": 5155.644931319356,
              "A_ds_shell": 138.2280912578236,
              "B_bot": 201.48451856399998,
              "B_mod_bot": 55.86894898700001,
              "h_drain": 28.665109974000018,
              "d_ds_perp": 3.510688710026443,
              "d_us_min": 21.101340232000013,
              "X_c_sec": 29.066374193,
              "Y_c_sec": 40.793548714999986,
              "W_sec_mean": 28.4627692655,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_core_eccentricity": 0.14426107971053515,
              "deep_soft_drainage_deficit": 0.18859251927106604
            },
            "loo_prediction_m": 0.15068445182380227,
            "observed_m": 0.1438072628529516,
            "observed_minus_loo_m": -0.006877188970850684
          },
          "matched_low_residual": {
            "sample_index": 20,
            "design_id": "lhs_020",
            "original_variables": {
              "upstream_slope": 1.4420337151416434,
              "downstream_slope": 1.6952251044089008,
              "crest_width": 10.899606022680448,
              "modulus_zone_height": 13.154257709251336,
              "secondary_top_width": 23.49373414724033,
              "secondary_height": 32.683541574470254,
              "secondary_bottom_width": 52.227206782238454,
              "secondary_right_slope_distance": 7.27158587734832
            },
            "geometric_primitives": {
              "A_tot": 8659.800493040262,
              "A_mod": 386.1615117911777,
              "A_sec": 1237.4142604653896,
              "A_drain": 1956.4091925520615,
              "A_main": 5079.815528231634,
              "A_ds_shell": 461.85458236732734,
              "B_bot": 235.65675416000002,
              "B_mod_bot": 50.418226704000006,
              "h_drain": 25.162200715999973,
              "d_ds_perp": 7.27158587734832,
              "d_us_min": 12.972606106,
              "X_c_sec": 29.586810307000007,
              "Y_c_sec": 39.436926298,
              "W_sec_mean": 37.8604704645,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_core_eccentricity": 0.1255504448088593,
              "deep_soft_drainage_deficit": 0.18096094215668337
            },
            "loo_prediction_m": 0.15141769150132883,
            "observed_m": 0.15120782641769512,
            "observed_minus_loo_m": -0.0002098650836337146
          }
        },
        {
          "high_residual": {
            "sample_index": 30,
            "design_id": "eval_0030",
            "original_variables": {
              "upstream_slope": 1.6313734091818333,
              "downstream_slope": 1.6791272994130848,
              "crest_width": 6.071929011493921,
              "modulus_zone_height": 19.592603002674878,
              "secondary_top_width": 14.385181814432144,
              "secondary_height": 20.009666811674833,
              "secondary_bottom_width": 61.474690763279796,
              "secondary_right_slope_distance": 6.464940689504147
            },
            "geometric_primitives": {
              "A_tot": 8790.413472151795,
              "A_mod": 725.9896508849653,
              "A_sec": 758.9653872971033,
              "A_drain": 3492.8760974358675,
              "A_main": 3812.5823365338565,
              "A_ds_shell": 250.80973158676687,
              "B_bot": 243.64897693799998,
              "B_mod_bot": 70.12193760200003,
              "h_drain": 31.39773018599999,
              "d_ds_perp": 6.464940689504147,
              "d_us_min": 43.306428987000004,
              "X_c_sec": 26.825576980999983,
              "Y_c_sec": 39.33241909899999,
              "W_sec_mean": 37.929936289000025,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_core_eccentricity": 0.11009928019449941,
              "deep_soft_drainage_deficit": 0.08362511409136761
            },
            "loo_prediction_m": 0.1364920972462234,
            "observed_m": 0.1425096541879413,
            "observed_minus_loo_m": 0.006017556941717878
          },
          "matched_low_residual": {
            "sample_index": 6,
            "design_id": "lhs_006",
            "original_variables": {
              "upstream_slope": 1.4984125475442482,
              "downstream_slope": 1.4983500043069138,
              "crest_width": 8.764640387706127,
              "modulus_zone_height": 15.56011237959704,
              "secondary_top_width": 18.525734632902335,
              "secondary_height": 25.22773911812784,
              "secondary_bottom_width": 42.48636761761311,
              "secondary_right_slope_distance": 5.681667919744493
            },
            "geometric_primitives": {
              "A_tot": 8124.688657783201,
              "A_mod": 459.94107055730245,
              "A_sec": 769.5986993184081,
              "A_drain": 2031.492431420449,
              "A_main": 4863.656456487041,
              "A_ds_shell": 254.52176138450596,
              "B_bot": 223.12506471600003,
              "B_mod_bot": 53.37986146999998,
              "h_drain": 30.21214850199999,
              "d_ds_perp": 5.681667919744493,
              "d_us_min": 24.910982525999998,
              "X_c_sec": 26.94885141000003,
              "Y_c_sec": 41.17478118700001,
              "W_sec_mean": 30.506051125,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_core_eccentricity": 0.12077913095196768,
              "deep_soft_drainage_deficit": 0.12187454201233337
            },
            "loo_prediction_m": 0.14123527375830952,
            "observed_m": 0.14069650531293582,
            "observed_minus_loo_m": -0.0005387684453737052
          }
        },
        {
          "high_residual": {
            "sample_index": 21,
            "design_id": "eval_0021",
            "original_variables": {
              "upstream_slope": 1.3067649867385627,
              "downstream_slope": 1.4604801163077354,
              "crest_width": 6.067888509482145,
              "modulus_zone_height": 15.058016548864543,
              "secondary_top_width": 14.974066410213709,
              "secondary_height": 22.15497637167573,
              "secondary_bottom_width": 47.734908745624125,
              "secondary_right_slope_distance": 9.691705060191453
            },
            "geometric_primitives": {
              "A_tot": 7306.629670651106,
              "A_mod": 356.5283997283624,
              "A_sec": 694.6579314617343,
              "A_drain": 3042.3260971476066,
              "A_main": 3213.1172423134026,
              "A_ds_shell": 376.6173875705929,
              "B_bot": 203.44402162600002,
              "B_mod_bot": 45.001169473999965,
              "h_drain": 33.78700707899998,
              "d_ds_perp": 9.691705060191453,
              "d_us_min": 13.165231336999994,
              "X_c_sec": 13.795161903000007,
              "Y_c_sec": 42.93543676899998,
              "W_sec_mean": 31.354487578000004,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_core_eccentricity": 0.06780814591032934,
              "deep_soft_drainage_deficit": 0.07803389574935964
            },
            "loo_prediction_m": 0.16246702560021098,
            "observed_m": 0.15663610636298597,
            "observed_minus_loo_m": -0.005830919237225007
          },
          "matched_low_residual": {
            "sample_index": 18,
            "design_id": "lhs_018",
            "original_variables": {
              "upstream_slope": 1.4073083504793016,
              "downstream_slope": 1.6098931341111617,
              "crest_width": 6.968367593672207,
              "modulus_zone_height": 14.064017121294864,
              "secondary_top_width": 25.924104461028058,
              "secondary_height": 23.459001636342418,
              "secondary_bottom_width": 33.408226510869504,
              "secondary_right_slope_distance": 9.609280043438371
            },
            "geometric_primitives": {
              "A_tot": 8053.077537985863,
              "A_mod": 361.127237449573,
              "A_sec": 695.9386246940131,
              "A_drain": 1913.1590371105046,
              "A_main": 5082.852638731771,
              "A_ds_shell": 423.96854022055777,
              "B_bot": 222.84127294400002,
              "B_mod_bot": 47.35163473,
              "h_drain": 33.47698124200002,
              "d_ds_perp": 9.609280043438371,
              "d_us_min": 3.4910248179999996,
              "X_c_sec": 17.13038473100002,
              "Y_c_sec": 44.713299483000014,
              "W_sec_mean": 29.666165486000025,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_core_eccentricity": 0.07687258515752997,
              "deep_soft_drainage_deficit": 0.1055648491461753
            },
            "loo_prediction_m": 0.14272627108090363,
            "observed_m": 0.14262663941530093,
            "observed_minus_loo_m": -9.963166560270076e-05
          }
        }
      ]
    },
    "dam_max_horizontal_displacement_m": {
      "current_formulas": [
        {
          "name": "soft_eccentric_area_fraction",
          "expression": "(X_c_sec / B_bot) * (A_sec / A_tot)",
          "physical_rationale": "将次堆石的相对面积与有符号偏心相乘以刻画软区对断面不对称的贡献。"
        },
        {
          "name": "main_drain_load_balance",
          "expression": "(A_main - A_drain) / A_tot",
          "physical_rationale": "候选假设：负残差高估子组中，高残断面相对对照表现为高刚度主堆石骨架面积相对底部排水支承更占优，即更大的主骨架与排水面积差占整体断面比例。主骨架占优可能把水平推力更多传给刚性承载路径，使实际水平位移小于模型预测；正残差和另一负残子组呈相反或较弱趋势，所以该解释仅作为特定刚度分配子组的候选假设。"
        },
        {
          "name": "cap_shell_over_soft_depth_factor",
          "expression": "((A_mod + A_ds_shell) / A_sec) * ((H_dam - Y_c_sec) / H_dam)",
          "physical_rationale": "候选假设：正残差两对中，高残差断面的 A_mod+A_ds_shell 相对 A_sec 及软区埋深因子更大；负残差中有一对帽盖和下游壳体相对软区明显减小，另有断面呈高位或顶部偏重软区。力学上，较大的顶部增模区与下游壳体在低位软区上方形成更宽的变形传递路径，可能增加向下游的相容位移；该解释只对深位软区子组作为候选假设成立。"
        }
      ],
      "loo_nrmse": 0.2791134664253034,
      "loo_formula_fallback": null,
      "top5_signs": {
        "underestimated": 3,
        "overestimated": 2,
        "zero": 0
      },
      "top5_direction": "mixed_or_tied",
      "matched_pairs": [
        {
          "high_residual": {
            "sample_index": 19,
            "design_id": "lhs_019",
            "original_variables": {
              "upstream_slope": 1.6745232365833786,
              "downstream_slope": 1.4008319126769921,
              "crest_width": 9.165674889626958,
              "modulus_zone_height": 16.760990775170942,
              "secondary_top_width": 24.924136868277312,
              "secondary_height": 36.21510900519744,
              "secondary_bottom_width": 57.11954572428753,
              "secondary_right_slope_distance": 10.172838528352186
            },
            "geometric_primitives": {
              "A_tot": 8367.75443834113,
              "A_mod": 547.9646538156717,
              "A_sec": 1485.6104541145198,
              "A_drain": 1491.4065154909613,
              "A_main": 4842.7728149199775,
              "A_ds_shell": 626.5786803350254,
              "B_bot": 229.341951425,
              "B_mod_bot": 59.01079680599997,
              "h_drain": 18.023900219999973,
              "d_ds_perp": 10.172838528352186,
              "d_us_min": 16.988913295000003,
              "X_c_sec": 21.663753582000027,
              "Y_c_sec": 33.76287860899998,
              "W_sec_mean": 41.02184129599998,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_eccentric_area_fraction": 0.016770508210026933,
              "main_drain_load_balance": 0.40050963781549614,
              "cap_shell_over_soft_depth_factor": 0.42989191871606736
            },
            "loo_prediction_m": 0.1070715246301364,
            "observed_m": 0.12216219021832611,
            "observed_minus_loo_m": 0.015090665588189708
          },
          "matched_low_residual": {
            "sample_index": 6,
            "design_id": "lhs_006",
            "original_variables": {
              "upstream_slope": 1.4984125475442482,
              "downstream_slope": 1.4983500043069138,
              "crest_width": 8.764640387706127,
              "modulus_zone_height": 15.56011237959704,
              "secondary_top_width": 18.525734632902335,
              "secondary_height": 25.22773911812784,
              "secondary_bottom_width": 42.48636761761311,
              "secondary_right_slope_distance": 5.681667919744493
            },
            "geometric_primitives": {
              "A_tot": 8124.688657783201,
              "A_mod": 459.94107055730245,
              "A_sec": 769.5986993184081,
              "A_drain": 2031.492431420449,
              "A_main": 4863.656456487041,
              "A_ds_shell": 254.52176138450596,
              "B_bot": 223.12506471600003,
              "B_mod_bot": 53.37986146999998,
              "h_drain": 30.21214850199999,
              "d_ds_perp": 5.681667919744493,
              "d_us_min": 24.910982525999998,
              "X_c_sec": 26.94885141000003,
              "Y_c_sec": 41.17478118700001,
              "W_sec_mean": 30.506051125,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_eccentric_area_fraction": 0.011440618342513025,
              "main_drain_load_balance": 0.3485873913892647,
              "cap_shell_over_soft_depth_factor": 0.41180463049638766
            },
            "loo_prediction_m": 0.06824999261219906,
            "observed_m": 0.07357472797703504,
            "observed_minus_loo_m": 0.005324735364835978
          }
        },
        {
          "high_residual": {
            "sample_index": 25,
            "design_id": "eval_0025",
            "original_variables": {
              "upstream_slope": 1.6882901702076196,
              "downstream_slope": 1.6950646217912435,
              "crest_width": 11.977469097822905,
              "modulus_zone_height": 18.36188363377005,
              "secondary_top_width": 11.75962993875146,
              "secondary_height": 27.23224811255932,
              "secondary_bottom_width": 25.882059396244586,
              "secondary_right_slope_distance": 13.828325916081667
            },
            "geometric_primitives": {
              "A_tot": 9408.289879923997,
              "A_mod": 766.9773790584494,
              "A_sec": 512.5339116636007,
              "A_drain": 1520.5571967507958,
              "A_main": 6608.221392451151,
              "A_ds_shell": 737.2519028341619,
              "B_bot": 254.94571919599997,
              "B_mod_bot": 73.42945398900002,
              "h_drain": 25.405868253999984,
              "d_ds_perp": 13.828325916081667,
              "d_us_min": 34.727601421,
              "X_c_sec": 30.76631473499998,
              "Y_c_sec": 37.31915736899998,
              "W_sec_mean": 18.820844667500026,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_eccentric_area_fraction": 0.006574150923512219,
              "main_drain_load_balance": 0.5407639709908103,
              "cap_shell_over_soft_depth_factor": 1.4547856780239483
            },
            "loo_prediction_m": 0.038527763397315756,
            "observed_m": 0.048193845100040454,
            "observed_minus_loo_m": 0.009666081702724698
          },
          "matched_low_residual": {
            "sample_index": 32,
            "design_id": "eval_0032",
            "original_variables": {
              "upstream_slope": 1.3320104010403158,
              "downstream_slope": 1.3892031960189342,
              "crest_width": 7.83259486220777,
              "modulus_zone_height": 17.462718952447176,
              "secondary_top_width": 18.06828137487173,
              "secondary_height": 24.233660977333784,
              "secondary_bottom_width": 33.185473438352346,
              "secondary_right_slope_distance": 10.66262403037399
            },
            "geometric_primitives": {
              "A_tot": 7306.114505716314,
              "A_mod": 494.31824134870135,
              "A_sec": 621.0330589964988,
              "A_drain": 1779.5526228276067,
              "A_main": 4411.210582543508,
              "A_ds_shell": 438.5025185213277,
              "B_bot": 201.80239653600003,
              "B_mod_bot": 52.634710216,
              "h_drain": 29.303620070000022,
              "d_ds_perp": 10.66262403037399,
              "d_us_min": 16.62668582699999,
              "X_c_sec": 19.091129449999983,
              "Y_c_sec": 40.229172285000004,
              "W_sec_mean": 25.6268774065,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_eccentric_area_fraction": 0.00804143485127701,
              "main_drain_load_balance": 0.360199386097342,
              "cap_shell_over_soft_depth_factor": 0.6854779179364483
            },
            "loo_prediction_m": 0.07587555668433751,
            "observed_m": 0.07496239422793384,
            "observed_minus_loo_m": -0.0009131624564036739
          }
        },
        {
          "high_residual": {
            "sample_index": 4,
            "design_id": "lhs_004",
            "original_variables": {
              "upstream_slope": 1.5014740020205382,
              "downstream_slope": 1.5240172039007247,
              "crest_width": 8.143891063356346,
              "modulus_zone_height": 15.05385059674546,
              "secondary_top_width": 20.880637571229904,
              "secondary_height": 37.75439049568013,
              "secondary_bottom_width": 25.12606122906331,
              "secondary_right_slope_distance": 11.414546318035475
            },
            "geometric_primitives": {
              "A_tot": 8159.104536977323,
              "A_mod": 428.5172402523009,
              "A_sec": 868.4774359454105,
              "A_drain": 856.1050285227655,
              "A_main": 6006.004832256849,
              "A_ds_shell": 777.3720588287769,
              "B_bot": 224.63023579199998,
              "B_mod_bot": 51.727709488000016,
              "h_drain": 18.191758907999997,
              "d_ds_perp": 11.414546318035475,
              "d_us_min": 10.47389064699999,
              "X_c_sec": 28.690286130999993,
              "Y_c_sec": 36.488301820000004,
              "W_sec_mean": 23.003349400000005,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_eccentric_area_fraction": 0.01359510945656136,
              "main_drain_load_balance": 0.6311844370169021,
              "cap_shell_over_soft_depth_factor": 0.7038561465467201
            },
            "loo_prediction_m": 0.09112213312268261,
            "observed_m": 0.08146443749533204,
            "observed_minus_loo_m": -0.009657695627350568
          },
          "matched_low_residual": {
            "sample_index": 14,
            "design_id": "lhs_014",
            "original_variables": {
              "upstream_slope": 1.4751002650383611,
              "downstream_slope": 1.595669529967744,
              "crest_width": 8.064570860520767,
              "modulus_zone_height": 19.724240431447242,
              "secondary_top_width": 28.442154455937654,
              "secondary_height": 33.84352919365724,
              "secondary_bottom_width": 38.69247372526842,
              "secondary_right_slope_distance": 14.13000948719288
            },
            "geometric_primitives": {
              "A_tot": 8277.12235269257,
              "A_mod": 713.7342492409941,
              "A_sec": 1136.036374381157,
              "A_drain": 1186.1341314208164,
              "A_main": 5241.217597649602,
              "A_ds_shell": 893.9669002097403,
              "B_bot": 227.901531182,
              "B_mod_bot": 67.11122239899998,
              "h_drain": 17.432230375000017,
              "d_ds_perp": 14.13000948719288,
              "d_us_min": 12.44616542999998,
              "X_c_sec": 24.94246610099998,
              "Y_c_sec": 33.492772806000005,
              "W_sec_mean": 33.567314090500005,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_eccentric_area_fraction": 0.015021213183172082,
              "main_drain_load_balance": 0.48991464586839845,
              "cap_shell_over_soft_depth_factor": 0.7746649047401243
            },
            "loo_prediction_m": 0.07348635566364638,
            "observed_m": 0.07637617202513365,
            "observed_minus_loo_m": 0.002889816361487274
          }
        },
        {
          "high_residual": {
            "sample_index": 1,
            "design_id": "lhs_001",
            "original_variables": {
              "upstream_slope": 1.6094092387020262,
              "downstream_slope": 1.398042551563892,
              "crest_width": 10.294596895348278,
              "modulus_zone_height": 17.22436933287436,
              "secondary_top_width": 22.75627311159097,
              "secondary_height": 28.75351770691302,
              "secondary_bottom_width": 46.6598200010875,
              "secondary_right_slope_distance": 12.621874746349105
            },
            "geometric_primitives": {
              "A_tot": 8262.411076914119,
              "A_mod": 581.5095874314916,
              "A_sec": 997.9784312431392,
              "A_drain": 1974.2684502131478,
              "A_main": 4708.654608026342,
              "A_ds_shell": 619.0967573010715,
              "B_bot": 225.44602486499997,
              "B_mod_bot": 60.22160053099998,
              "h_drain": 25.022112960000015,
              "d_ds_perp": 12.621874746349105,
              "d_us_min": 16.09770578199999,
              "X_c_sec": 15.802945496999996,
              "Y_c_sec": 37.74865155499998,
              "W_sec_mean": 34.70804655650002,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_eccentric_area_fraction": 0.008466615102115866,
              "main_drain_load_balance": 0.3309428848745256,
              "cap_shell_over_soft_depth_factor": 0.589348150350608
            },
            "loo_prediction_m": 0.08975770735048216,
            "observed_m": 0.08021133624444485,
            "observed_minus_loo_m": -0.009546371106037316
          },
          "matched_low_residual": {
            "sample_index": 18,
            "design_id": "lhs_018",
            "original_variables": {
              "upstream_slope": 1.4073083504793016,
              "downstream_slope": 1.6098931341111617,
              "crest_width": 6.968367593672207,
              "modulus_zone_height": 14.064017121294864,
              "secondary_top_width": 25.924104461028058,
              "secondary_height": 23.459001636342418,
              "secondary_bottom_width": 33.408226510869504,
              "secondary_right_slope_distance": 9.609280043438371
            },
            "geometric_primitives": {
              "A_tot": 8053.077537985863,
              "A_mod": 361.127237449573,
              "A_sec": 695.9386246940131,
              "A_drain": 1913.1590371105046,
              "A_main": 5082.852638731771,
              "A_ds_shell": 423.96854022055777,
              "B_bot": 222.84127294400002,
              "B_mod_bot": 47.35163473,
              "h_drain": 33.47698124200002,
              "d_ds_perp": 9.609280043438371,
              "d_us_min": 3.4910248179999996,
              "X_c_sec": 17.13038473100002,
              "Y_c_sec": 44.713299483000014,
              "W_sec_mean": 29.666165486000025,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_eccentric_area_fraction": 0.0066432492347000575,
              "main_drain_load_balance": 0.3936002834531295,
              "cap_shell_over_soft_depth_factor": 0.44646809324243153
            },
            "loo_prediction_m": 0.0532053588997993,
            "observed_m": 0.049864421905795264,
            "observed_minus_loo_m": -0.0033409369940040345
          }
        },
        {
          "high_residual": {
            "sample_index": 11,
            "design_id": "lhs_011",
            "original_variables": {
              "upstream_slope": 1.3717017078040479,
              "downstream_slope": 1.320912292949585,
              "crest_width": 9.763338568555401,
              "modulus_zone_height": 5.227712698589564,
              "secondary_top_width": 14.210190616724496,
              "secondary_height": 29.327369356524045,
              "secondary_bottom_width": 49.60905740491252,
              "secondary_right_slope_distance": 2.713831948671564
            },
            "geometric_primitives": {
              "A_tot": 7364.980022583283,
              "A_mod": 68.12974081013377,
              "A_sec": 935.825329374874,
              "A_drain": 2780.5265482474133,
              "A_main": 3580.498404150864,
              "A_ds_shell": 126.00460934764087,
              "B_bot": 201.616770115,
              "B_mod_bot": 20.240492020000005,
              "h_drain": 36.444917944999986,
              "d_ds_perp": 2.713831948671564,
              "d_us_min": 1.906027913999992,
              "X_c_sec": 16.229618721999998,
              "Y_c_sec": 48.39741523800001,
              "W_sec_mean": 31.909624011000005,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_eccentric_area_fraction": 0.010228333658466019,
              "main_drain_load_balance": 0.10861833344428529,
              "cap_shell_over_soft_depth_factor": 0.07177276598346778
            },
            "loo_prediction_m": 0.08757528382680994,
            "observed_m": 0.09700693589355729,
            "observed_minus_loo_m": 0.009431652066747342
          },
          "matched_low_residual": {
            "sample_index": 10,
            "design_id": "lhs_010",
            "original_variables": {
              "upstream_slope": 1.5468086093385431,
              "downstream_slope": 1.3012624657189515,
              "crest_width": 11.291715983318776,
              "modulus_zone_height": 13.534298401099962,
              "secondary_top_width": 26.23483907948985,
              "secondary_height": 27.391521289178613,
              "secondary_bottom_width": 63.350919984777626,
              "secondary_right_slope_distance": 8.592072960145197
            },
            "geometric_primitives": {
              "A_tot": 7897.915705416024,
              "A_mod": 372.63251196399375,
              "A_sec": 1226.9451132903607,
              "A_drain": 2947.1016370375714,
              "A_main": 3351.236443124097,
              "A_ds_shell": 381.73199744837643,
              "B_bot": 214.648971028,
              "B_mod_bot": 47.245817955999996,
              "h_drain": 30.074180309999974,
              "d_ds_perp": 8.592072960145197,
              "d_us_min": 7.236450769000015,
              "X_c_sec": 9.331825449999997,
              "Y_c_sec": 41.878520691000006,
              "W_sec_mean": 44.792879532000015,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_eccentric_area_fraction": 0.006753834853107899,
              "main_drain_load_balance": 0.0511698049409908,
              "cap_shell_over_soft_depth_factor": 0.26688239778945955
            },
            "loo_prediction_m": 0.09524170200406523,
            "observed_m": 0.09705632591519942,
            "observed_minus_loo_m": 0.0018146239111341833
          }
        }
      ]
    }
  }
}
```
