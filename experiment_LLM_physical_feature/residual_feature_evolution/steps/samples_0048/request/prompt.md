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
  "sample_count": 48,
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
          "name": "elevated_drain_soft_cap_deficit",
          "expression": "(A_drain / A_tot) * (Y_c_sec / H_dam) * (A_sec / (A_mod + A_sec))",
          "physical_rationale": "候选假设：高绝对残差中的多数负残差断面以及一个正残差断面，相对对照表现为底部排水面积占比更大、软区形心更高、顶部帽盖相对软区更弱，因此该指标增大；其他正残差断面不完全支持。力学上，高位软区上方帽盖不足而底部排水刚体支承占比大，会改变竖向拱架与压缩柱路径：在排水支承占主导时模型可能把坝顶沉降算得过高，在帽盖极薄时又可能低估局部沉降，所以此解释仅作为符号混合子组的候选假设。"
        }
      ],
      "loo_nrmse": 0.21119342737644065,
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
              "narrow_crest_broad_base_index": 25.02182918189172,
              "soft_area_ratio": 0.17753992006594252,
              "elevated_drain_soft_cap_deficit": 0.059407268784753164
            },
            "loo_prediction_m": 0.07633653010107609,
            "observed_m": 0.08383617259362591,
            "observed_minus_loo_m": 0.007499642492549821
          },
          "matched_low_residual": {
            "sample_index": 40,
            "design_id": "eval_0040",
            "original_variables": {
              "upstream_slope": 1.663573120906949,
              "downstream_slope": 1.674020792171359,
              "crest_width": 10.857982898131013,
              "modulus_zone_height": 5.0292914267629385,
              "secondary_top_width": 16.074523124843836,
              "secondary_height": 37.080765683203936,
              "secondary_bottom_width": 67.6366621768102,
              "secondary_right_slope_distance": 2.9419813733547926
            },
            "geometric_primitives": {
              "A_tot": 9203.913581336039,
              "A_mod": 87.56153094098183,
              "A_sec": 1552.0374236138073,
              "A_drain": 2705.388622663073,
              "A_main": 4858.926004118179,
              "A_ds_shell": 205.7830245656496,
              "B_bot": 250.439927957,
              "B_mod_bot": 25.966692556000027,
              "h_drain": 28.889942889999986,
              "d_ds_perp": 2.9419813733547926,
              "d_us_min": 4.543058373000008,
              "X_c_sec": 26.59827887199998,
              "Y_c_sec": 43.623658818000024,
              "W_sec_mean": 41.855592650499986,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "narrow_crest_broad_base_index": 23.06505087608015,
              "soft_area_ratio": 0.16862798742059829,
              "elevated_drain_soft_cap_deficit": 0.16402575732749256
            },
            "loo_prediction_m": 0.09175174712358966,
            "observed_m": 0.09257967335395396,
            "observed_minus_loo_m": 0.0008279262303643037
          }
        },
        {
          "high_residual": {
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
              "narrow_crest_broad_base_index": 32.43139618336402,
              "soft_area_ratio": 0.08109905499432402,
              "elevated_drain_soft_cap_deficit": 0.09998456837692517
            },
            "loo_prediction_m": 0.0752087702475619,
            "observed_m": 0.08108230883414239,
            "observed_minus_loo_m": 0.00587353858658049
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
              "narrow_crest_broad_base_index": 31.978977852195452,
              "soft_area_ratio": 0.08641896485055733,
              "elevated_drain_soft_cap_deficit": 0.09450680576533416
            },
            "loo_prediction_m": 0.07486677570897991,
            "observed_m": 0.07408918199300787,
            "observed_minus_loo_m": -0.0007775937159720386
          }
        },
        {
          "high_residual": {
            "sample_index": 45,
            "design_id": "eval_0045",
            "original_variables": {
              "upstream_slope": 1.6643347464501859,
              "downstream_slope": 1.4260760571807622,
              "crest_width": 11.799292227253318,
              "modulus_zone_height": 18.962352043017745,
              "secondary_top_width": 12.44045464321971,
              "secondary_height": 20.33594235777855,
              "secondary_bottom_width": 42.75813332758844,
              "secondary_right_slope_distance": 3.3349734395742416
            },
            "geometric_primitives": {
              "A_tot": 8595.69245209823,
              "A_mod": 738.9806991282018,
              "A_sec": 561.2576516042626,
              "A_drain": 2294.843520282626,
              "A_main": 5000.6105810831405,
              "A_ds_shell": 115.94895290367178,
              "B_bot": 233.089687186,
              "B_mod_bot": 68.888180605,
              "h_drain": 31.70170559899998,
              "d_ds_perp": 3.3349734395742416,
              "d_us_min": 50.85590718899999,
              "X_c_sec": 32.99781848499998,
              "Y_c_sec": 40.008099671000025,
              "W_sec_mean": 27.599293985500026,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "narrow_crest_broad_base_index": 19.754548213292235,
              "soft_area_ratio": 0.06529522254687674,
              "elevated_drain_soft_cap_deficit": 0.06230569088609921
            },
            "loo_prediction_m": 0.057376379233568975,
            "observed_m": 0.06226105178888307,
            "observed_minus_loo_m": 0.004884672555314093
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
              "elevated_drain_soft_cap_deficit": 0.03900158886449066
            },
            "loo_prediction_m": 0.07381538141666703,
            "observed_m": 0.07307196200591454,
            "observed_minus_loo_m": -0.0007434194107524966
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
              "elevated_drain_soft_cap_deficit": 0.2158514601795583
            },
            "loo_prediction_m": 0.0842423763484217,
            "observed_m": 0.07951806597548103,
            "observed_minus_loo_m": -0.004724310372940674
          },
          "matched_low_residual": {
            "sample_index": 41,
            "design_id": "eval_0041",
            "original_variables": {
              "upstream_slope": 1.3161638882011175,
              "downstream_slope": 1.4899378784000874,
              "crest_width": 6.085398696362972,
              "modulus_zone_height": 14.327679295092821,
              "secondary_top_width": 25.771700870245695,
              "secondary_height": 24.885670598596334,
              "secondary_bottom_width": 46.73315422143787,
              "secondary_right_slope_distance": 5.833214187994599
            },
            "geometric_primitives": {
              "A_tot": 7413.979788183631,
              "A_mod": 330.327076975829,
              "A_sec": 902.1659703276371,
              "A_drain": 2249.0214665753465,
              "A_main": 3932.4652743048196,
              "A_ds_shell": 256.4097486659593,
              "B_bot": 206.33692491600002,
              "B_mod_bot": 43.62346387500003,
              "h_drain": 31.786650106000025,
              "d_ds_perp": 5.833214187994599,
              "d_us_min": 7.701682787999999,
              "X_c_sec": 19.753220667999983,
              "Y_c_sec": 43.03039369700002,
              "W_sec_mean": 36.25242754599999,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "narrow_crest_broad_base_index": 33.90688683049154,
              "soft_area_ratio": 0.12168443887121264,
              "elevated_drain_soft_cap_deficit": 0.12911829016462514
            },
            "loo_prediction_m": 0.08059067306455038,
            "observed_m": 0.08107417861630467,
            "observed_minus_loo_m": 0.0004835055517542902
          }
        },
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
              "elevated_drain_soft_cap_deficit": 0.1560034604301601
            },
            "loo_prediction_m": 0.08167017404171102,
            "observed_m": 0.07708631363929866,
            "observed_minus_loo_m": -0.0045838604024123625
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
              "elevated_drain_soft_cap_deficit": 0.07372873442557262
            },
            "loo_prediction_m": 0.0742008012914495,
            "observed_m": 0.07437356009175546,
            "observed_minus_loo_m": 0.00017275880030595547
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
      "loo_nrmse": 0.25392301786884786,
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
            "sample_index": 31,
            "design_id": "eval_0031",
            "original_variables": {
              "upstream_slope": 1.3630816973745823,
              "downstream_slope": 1.3697291061282157,
              "crest_width": 6.585983991622925,
              "modulus_zone_height": 14.717251099646091,
              "secondary_top_width": 12.489079479128122,
              "secondary_height": 25.15388460829854,
              "secondary_bottom_width": 41.36039569508284,
              "secondary_right_slope_distance": 6.0032567055895925
            },
            "geometric_primitives": {
              "A_tot": 7249.350975337761,
              "A_mod": 343.7334144502847,
              "A_sec": 677.2617423696712,
              "A_drain": 2165.501912655367,
              "A_main": 4062.8539058624388,
              "A_ds_shell": 251.86419530972142,
              "B_bot": 201.41397894099998,
              "B_mod_bot": 43.944044393999974,
              "h_drain": 31.128864292000003,
              "d_ds_perp": 6.0032567055895925,
              "d_us_min": 21.600102883999995,
              "X_c_sec": 22.612904904999994,
              "Y_c_sec": 41.45810361100001,
              "W_sec_mean": 26.92473758699998,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_core_eccentricity": 0.11227078191838895,
              "deep_soft_drainage_deficit": 0.10476758365273954
            },
            "loo_prediction_m": 0.15813600923174106,
            "observed_m": 0.14523452720601138,
            "observed_minus_loo_m": -0.012901482025729688
          },
          "matched_low_residual": {
            "sample_index": 38,
            "design_id": "eval_0038",
            "original_variables": {
              "upstream_slope": 1.3465939488261938,
              "downstream_slope": 1.4056917678564786,
              "crest_width": 6.198216548189521,
              "modulus_zone_height": 10.847343560308218,
              "secondary_top_width": 16.56683437526226,
              "secondary_height": 24.447562992572784,
              "secondary_bottom_width": 34.648930304683745,
              "secondary_right_slope_distance": 3.9028568333014846
            },
            "geometric_primitives": {
              "A_tot": 7275.02454225306,
              "A_mod": 192.1977331747983,
              "A_sec": 626.0503166214861,
              "A_drain": 1955.5371611030084,
              "A_main": 4501.239331353769,
              "A_ds_shell": 160.2917079005102,
              "B_bot": 202.467355073,
              "B_mod_bot": 32.99859102000002,
              "h_drain": 35.705093447000024,
              "d_ds_perp": 3.9028568333014846,
              "d_us_min": 10.02189313400001,
              "X_c_sec": 21.332502328000032,
              "Y_c_sec": 46.490310137999984,
              "W_sec_mean": 25.607882340000003,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_core_eccentricity": 0.10536267597464567,
              "deep_soft_drainage_deficit": 0.09015220807786412
            },
            "loo_prediction_m": 0.1508127406356134,
            "observed_m": 0.1518017469838714,
            "observed_minus_loo_m": 0.0009890063482579925
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
              "soft_core_eccentricity": 0.12067790285722395,
              "deep_soft_drainage_deficit": 0.12496067101678475
            },
            "loo_prediction_m": 0.1210533607600759,
            "observed_m": 0.1320332795531413,
            "observed_minus_loo_m": 0.010979918793065382
          },
          "matched_low_residual": {
            "sample_index": 46,
            "design_id": "eval_0046",
            "original_variables": {
              "upstream_slope": 1.3095830541104079,
              "downstream_slope": 1.4183809250593185,
              "crest_width": 7.634220125153661,
              "modulus_zone_height": 15.260290312580764,
              "secondary_top_width": 27.250584866851568,
              "secondary_height": 22.001241110265255,
              "secondary_bottom_width": 31.37393809389323,
              "secondary_right_slope_distance": 3.7364539681002498
            },
            "geometric_primitives": {
              "A_tot": 7310.475739232889,
              "A_mod": 383.2194229221109,
              "A_sec": 644.9061323080062,
              "A_drain": 1389.8307554752264,
              "A_main": 4892.519428527544,
              "A_ds_shell": 139.27214511441463,
              "B_bot": 202.10355007400003,
              "B_mod_bot": 46.423110080999976,
              "h_drain": 33.73846857699999,
              "d_ds_perp": 3.7364539681002498,
              "d_us_min": 12.978903067000005,
              "X_c_sec": 24.171861794999984,
              "Y_c_sec": 44.48117974500002,
              "W_sec_mean": 29.312261480499984,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_core_eccentricity": 0.1196013716045536,
              "deep_soft_drainage_deficit": 0.12643156824504517
            },
            "loo_prediction_m": 0.1449459719971893,
            "observed_m": 0.14458871281002406,
            "observed_minus_loo_m": -0.0003572591871652542
          }
        },
        {
          "high_residual": {
            "sample_index": 9,
            "design_id": "lhs_009",
            "original_variables": {
              "upstream_slope": 1.521905918548933,
              "downstream_slope": 1.4341577509152494,
              "crest_width": 6.2300210381399586,
              "modulus_zone_height": 9.880374085179792,
              "secondary_top_width": 11.206422680062811,
              "secondary_height": 35.685761868782116,
              "secondary_bottom_width": 62.945727953476855,
              "secondary_right_slope_distance": 6.127392062091551
            },
            "geometric_primitives": {
              "A_tot": 7833.748071896116,
              "A_mod": 177.90726934915705,
              "A_sec": 1323.0879947750243,
              "A_drain": 2345.947228942543,
              "A_main": 3986.80557882939,
              "A_ds_shell": 374.8224797625352,
              "B_bot": 217.57872806899996,
              "B_mod_bot": 32.93087933900003,
              "h_drain": 25.433864045999997,
              "d_ds_perp": 6.127392062091551,
              "d_us_min": 11.431447804000015,
              "X_c_sec": 20.587154639999994,
              "Y_c_sec": 39.12681695700002,
              "W_sec_mean": 37.076075316499995,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_core_eccentricity": 0.0946193353675239,
              "deep_soft_drainage_deficit": 0.16994043175725815
            },
            "loo_prediction_m": 0.1763994873201386,
            "observed_m": 0.1658134202298651,
            "observed_minus_loo_m": -0.010586067090273515
          },
          "matched_low_residual": {
            "sample_index": 40,
            "design_id": "eval_0040",
            "original_variables": {
              "upstream_slope": 1.663573120906949,
              "downstream_slope": 1.674020792171359,
              "crest_width": 10.857982898131013,
              "modulus_zone_height": 5.0292914267629385,
              "secondary_top_width": 16.074523124843836,
              "secondary_height": 37.080765683203936,
              "secondary_bottom_width": 67.6366621768102,
              "secondary_right_slope_distance": 2.9419813733547926
            },
            "geometric_primitives": {
              "A_tot": 9203.913581336039,
              "A_mod": 87.56153094098183,
              "A_sec": 1552.0374236138073,
              "A_drain": 2705.388622663073,
              "A_main": 4858.926004118179,
              "A_ds_shell": 205.7830245656496,
              "B_bot": 250.439927957,
              "B_mod_bot": 25.966692556000027,
              "h_drain": 28.889942889999986,
              "d_ds_perp": 2.9419813733547926,
              "d_us_min": 4.543058373000008,
              "X_c_sec": 26.59827887199998,
              "Y_c_sec": 43.623658818000024,
              "W_sec_mean": 41.855592650499986,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_core_eccentricity": 0.10620622314093162,
              "deep_soft_drainage_deficit": 0.1496438312170615
            },
            "loo_prediction_m": 0.1549645842716876,
            "observed_m": 0.15472294791486899,
            "observed_minus_loo_m": -0.00024163635681861972
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
              "soft_core_eccentricity": 0.09446049206171757,
              "deep_soft_drainage_deficit": 0.2713431235989366
            },
            "loo_prediction_m": 0.18934477920747667,
            "observed_m": 0.19813934916032805,
            "observed_minus_loo_m": 0.008794569952851383
          },
          "matched_low_residual": {
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
              "soft_core_eccentricity": 0.12772228115170672,
              "deep_soft_drainage_deficit": 0.2552757608238603
            },
            "loo_prediction_m": 0.15936273008943042,
            "observed_m": 0.16155690587457133,
            "observed_minus_loo_m": 0.002194175785140917
          }
        },
        {
          "high_residual": {
            "sample_index": 41,
            "design_id": "eval_0041",
            "original_variables": {
              "upstream_slope": 1.3161638882011175,
              "downstream_slope": 1.4899378784000874,
              "crest_width": 6.085398696362972,
              "modulus_zone_height": 14.327679295092821,
              "secondary_top_width": 25.771700870245695,
              "secondary_height": 24.885670598596334,
              "secondary_bottom_width": 46.73315422143787,
              "secondary_right_slope_distance": 5.833214187994599
            },
            "geometric_primitives": {
              "A_tot": 7413.979788183631,
              "A_mod": 330.327076975829,
              "A_sec": 902.1659703276371,
              "A_drain": 2249.0214665753465,
              "A_main": 3932.4652743048196,
              "A_ds_shell": 256.4097486659593,
              "B_bot": 206.33692491600002,
              "B_mod_bot": 43.62346387500003,
              "h_drain": 31.786650106000025,
              "d_ds_perp": 5.833214187994599,
              "d_us_min": 7.701682787999999,
              "X_c_sec": 19.753220667999983,
              "Y_c_sec": 43.03039369700002,
              "W_sec_mean": 36.25242754599999,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_core_eccentricity": 0.09573284411426379,
              "deep_soft_drainage_deficit": 0.11981637330180651
            },
            "loo_prediction_m": 0.15523157314683783,
            "observed_m": 0.16332233183291922,
            "observed_minus_loo_m": 0.00809075868608139
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
              "soft_core_eccentricity": 0.09460308587858751,
              "deep_soft_drainage_deficit": 0.1180612815936371
            },
            "loo_prediction_m": 0.15709452049607645,
            "observed_m": 0.15799288529580513,
            "observed_minus_loo_m": 0.0008983647997286848
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
          "name": "upstream_clearance_asymmetry_index",
          "expression": "(upstream_slope - downstream_slope) * (d_us_min / H_dam)",
          "physical_rationale": "候选假设：五个高残差断面的该指标均高于对照；它把上下游坡比不对称与软区上游净距结合，而不是单纯使用总体底宽。正残差的 lhs_019、eval_0034 中，更大的不对称和有效上游净距可能使自重合力与水平推力线偏移，模型低估水平位移；负残差的 lhs_001、eval_0031、lhs_002 中，即使主骨架或排水占优，局部不对称仍改变剪切路径，模型可能对约束估计过强而高估。该式是上游约束—推力线偏移的候选机制。"
        },
        {
          "name": "drain_slope_asymmetry_shear",
          "expression": "abs(upstream_slope - downstream_slope) * (A_drain / A_tot)",
          "physical_rationale": "候选假设：正残差 lhs_019、eval_0037、eval_0034 的坡比不对称与底部排水面积占比组合均高于对照，负残差 lhs_017 亦高，负残差 lhs_001 不服从。该无量纲量表征坡比不对称引起的推力偏移与刚性底部排水约束的耦合：排水越强，偏移推力越可能沿底部剪切带转化为下游水平位移；正残差子组中模型可能低估该剪切转化，负残差子组中若主骨架很厚则可能高估。"
        }
      ],
      "loo_nrmse": 0.23482826338494311,
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
              "upstream_clearance_asymmetry_index": 0.06283402934377516,
              "drain_slope_asymmetry_shear": 0.04878071252151286
            },
            "loo_prediction_m": 0.10945739301907226,
            "observed_m": 0.12216219021832611,
            "observed_minus_loo_m": 0.012704797199253853
          },
          "matched_low_residual": {
            "sample_index": 9,
            "design_id": "lhs_009",
            "original_variables": {
              "upstream_slope": 1.521905918548933,
              "downstream_slope": 1.4341577509152494,
              "crest_width": 6.2300210381399586,
              "modulus_zone_height": 9.880374085179792,
              "secondary_top_width": 11.206422680062811,
              "secondary_height": 35.685761868782116,
              "secondary_bottom_width": 62.945727953476855,
              "secondary_right_slope_distance": 6.127392062091551
            },
            "geometric_primitives": {
              "A_tot": 7833.748071896116,
              "A_mod": 177.90726934915705,
              "A_sec": 1323.0879947750243,
              "A_drain": 2345.947228942543,
              "A_main": 3986.80557882939,
              "A_ds_shell": 374.8224797625352,
              "B_bot": 217.57872806899996,
              "B_mod_bot": 32.93087933900003,
              "h_drain": 25.433864045999997,
              "d_ds_perp": 6.127392062091551,
              "d_us_min": 11.431447804000015,
              "X_c_sec": 20.587154639999994,
              "Y_c_sec": 39.12681695700002,
              "W_sec_mean": 37.076075316499995,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_eccentric_area_fraction": 0.01598081857488956,
              "upstream_clearance_asymmetry_index": 0.013555251327041851,
              "drain_slope_asymmetry_shear": 0.026277660299484223
            },
            "loo_prediction_m": 0.08647646699272973,
            "observed_m": 0.08738836831204541,
            "observed_minus_loo_m": 0.0009119013193156777
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
              "upstream_clearance_asymmetry_index": 0.04597998299548281,
              "drain_slope_asymmetry_shear": 0.05050518280418756
            },
            "loo_prediction_m": 0.090430714497061,
            "observed_m": 0.08021133624444485,
            "observed_minus_loo_m": -0.010219378252616151
          },
          "matched_low_residual": {
            "sample_index": 16,
            "design_id": "lhs_016",
            "original_variables": {
              "upstream_slope": 1.5904733072542485,
              "downstream_slope": 1.3508043169679371,
              "crest_width": 11.693416459617298,
              "modulus_zone_height": 18.970891762402378,
              "secondary_top_width": 29.273667476470944,
              "secondary_height": 34.417047583299556,
              "secondary_bottom_width": 29.75174025290207,
              "secondary_right_slope_distance": 14.580711971134479
            },
            "geometric_primitives": {
              "A_tot": 8180.935829827038,
              "A_mod": 702.2368509688888,
              "A_sec": 1015.7401332483503,
              "A_drain": 959.5100318499856,
              "A_main": 5503.448813759815,
              "A_ds_shell": 836.5874225057545,
              "B_bot": 221.947956142,
              "B_mod_bot": 65.53266242199996,
              "h_drain": 17.612060654000004,
              "d_ds_perp": 14.580711971134479,
              "d_us_min": 12.149648116999998,
              "X_c_sec": 19.369528333999995,
              "Y_c_sec": 34.774124641000014,
              "W_sec_mean": 29.512703864500025,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_eccentric_area_fraction": 0.010835464198331911,
              "upstream_clearance_asymmetry_index": 0.0393499175207483,
              "drain_slope_asymmetry_shear": 0.028109840400488077
            },
            "loo_prediction_m": 0.11595421454719895,
            "observed_m": 0.11281189766268739,
            "observed_minus_loo_m": -0.0031423168845115673
          }
        },
        {
          "high_residual": {
            "sample_index": 17,
            "design_id": "lhs_017",
            "original_variables": {
              "upstream_slope": 1.6844315309886524,
              "downstream_slope": 1.4551043198867881,
              "crest_width": 6.360159107004971,
              "modulus_zone_height": 7.447526494976314,
              "secondary_top_width": 10.016459459486484,
              "secondary_height": 38.77411950955866,
              "secondary_bottom_width": 37.57497790629451,
              "secondary_right_slope_distance": 3.53078965948255
            },
            "geometric_primitives": {
              "A_tot": 8343.997175873314,
              "A_mod": 116.8886371461942,
              "A_sec": 922.65804003732,
              "A_drain": 1307.3043572125432,
              "A_main": 5997.146141477255,
              "A_ds_shell": 233.24576754904328,
              "B_bot": 231.285807562,
              "B_mod_bot": 27.62799936600001,
              "h_drain": 24.778353994999975,
              "d_ds_perp": 3.53078965948255,
              "d_us_min": 11.816037340999998,
              "X_c_sec": 32.28448626599999,
              "Y_c_sec": 40.423293504000014,
              "W_sec_mean": 23.795718683499985,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_eccentric_area_fraction": 0.015435173594292478,
              "upstream_clearance_asymmetry_index": 0.03661809310387861,
              "drain_slope_asymmetry_shear": 0.03593007715387795
            },
            "loo_prediction_m": 0.09049608228555081,
            "observed_m": 0.0807759804523365,
            "observed_minus_loo_m": -0.009720101833214312
          },
          "matched_low_residual": {
            "sample_index": 12,
            "design_id": "lhs_012",
            "original_variables": {
              "upstream_slope": 1.4358822591610336,
              "downstream_slope": 1.3705155660015915,
              "crest_width": 11.72827352663696,
              "modulus_zone_height": 10.858807129173226,
              "secondary_top_width": 13.057879130529466,
              "secondary_height": 39.41897416299703,
              "secondary_bottom_width": 43.37466410070071,
              "secondary_right_slope_distance": 10.671255276130882
            },
            "geometric_primitives": {
              "A_tot": 7815.085045741038,
              "A_mod": 257.57908015435964,
              "A_sec": 1112.2564817836796,
              "A_drain": 1441.4633632943508,
              "A_main": 5003.786120508648,
              "A_ds_shell": 704.5739071908122,
              "B_bot": 212.00170807900003,
              "B_mod_bot": 39.310873696000016,
              "h_drain": 20.722218708000014,
              "d_ds_perp": 10.671255276130882,
              "d_us_min": 8.612333860000007,
              "X_c_sec": 22.910333228000013,
              "Y_c_sec": 36.90225128499998,
              "W_sec_mean": 28.216271616,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_eccentric_area_fraction": 0.01538024603382966,
              "upstream_clearance_asymmetry_index": 0.007607564659639108,
              "drain_slope_asymmetry_shear": 0.012056643378486076
            },
            "loo_prediction_m": 0.10407220855223,
            "observed_m": 0.10412201094086328,
            "observed_minus_loo_m": 4.9802388633279904e-05
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
              "upstream_clearance_asymmetry_index": -0.0031791953302966565,
              "drain_slope_asymmetry_shear": 0.0010948792225779441
            },
            "loo_prediction_m": 0.03900303261669353,
            "observed_m": 0.048193845100040454,
            "observed_minus_loo_m": 0.009190812483346925
          },
          "matched_low_residual": {
            "sample_index": 47,
            "design_id": "eval_0047",
            "original_variables": {
              "upstream_slope": 1.3592205163091422,
              "downstream_slope": 1.386825466528535,
              "crest_width": 6.003654453903437,
              "modulus_zone_height": 14.168338072486222,
              "secondary_top_width": 27.092327158898115,
              "secondary_height": 21.997824627906084,
              "secondary_bottom_width": 25.615036729723215,
              "secondary_right_slope_distance": 3.7122870655730367
            },
            "geometric_primitives": {
              "A_tot": 7244.183433494947,
              "A_mod": 313.67380411541694,
              "A_sec": 579.7236737088316,
              "A_drain": 1078.822039547066,
              "A_main": 5271.963916123632,
              "A_ds_shell": 136.10749355327542,
              "B_bot": 201.81105267400002,
              "B_mod_bot": 42.05317295399999,
              "h_drain": 34.83383730000003,
              "d_ds_perp": 3.7122870655730367,
              "d_us_min": 8.908790076000003,
              "X_c_sec": 22.037229592000017,
              "Y_c_sec": 45.93550937200001,
              "W_sec_mean": 26.353681943999987,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_eccentric_area_fraction": 0.00873863591756307,
              "upstream_clearance_asymmetry_index": -0.0033233338724729684,
              "drain_slope_asymmetry_shear": 0.004110998702708543
            },
            "loo_prediction_m": 0.08118839979809482,
            "observed_m": 0.08426405409010362,
            "observed_minus_loo_m": 0.003075654292008803
          }
        },
        {
          "high_residual": {
            "sample_index": 45,
            "design_id": "eval_0045",
            "original_variables": {
              "upstream_slope": 1.6643347464501859,
              "downstream_slope": 1.4260760571807622,
              "crest_width": 11.799292227253318,
              "modulus_zone_height": 18.962352043017745,
              "secondary_top_width": 12.44045464321971,
              "secondary_height": 20.33594235777855,
              "secondary_bottom_width": 42.75813332758844,
              "secondary_right_slope_distance": 3.3349734395742416
            },
            "geometric_primitives": {
              "A_tot": 8595.69245209823,
              "A_mod": 738.9806991282018,
              "A_sec": 561.2576516042626,
              "A_drain": 2294.843520282626,
              "A_main": 5000.6105810831405,
              "A_ds_shell": 115.94895290367178,
              "B_bot": 233.089687186,
              "B_mod_bot": 68.888180605,
              "h_drain": 31.70170559899998,
              "d_ds_perp": 3.3349734395742416,
              "d_us_min": 50.85590718899999,
              "X_c_sec": 32.99781848499998,
              "Y_c_sec": 40.008099671000025,
              "W_sec_mean": 27.599293985500026,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_eccentric_area_fraction": 0.00924365177864003,
              "upstream_clearance_asymmetry_index": 0.16374137551971077,
              "drain_slope_asymmetry_shear": 0.06360934994684475
            },
            "loo_prediction_m": 0.07907696277488,
            "observed_m": 0.07041981318398705,
            "observed_minus_loo_m": -0.008657149590892949
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
              "upstream_clearance_asymmetry_index": 0.024011926754337592,
              "drain_slope_asymmetry_shear": 0.09162536912533971
            },
            "loo_prediction_m": 0.09950876785824846,
            "observed_m": 0.09705632591519942,
            "observed_minus_loo_m": -0.0024524419430490413
          }
        }
      ]
    }
  }
}
```
