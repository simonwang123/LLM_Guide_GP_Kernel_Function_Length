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
  "sample_count": 56,
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
          "name": "combined_cap_shell_area_fraction",
          "expression": "(A_mod + A_ds_shell) / A_tot",
          "physical_rationale": "候选假设：三条正残差低估断面的顶部增模帽盖加下游壳体面积相对坝体总面积均高于对照，尤其 lhs_019 与 eval_0028；两条负残差高估断面则低于对照。该量比单一帽盖面积更稳定，代表顶部帽盖与下游壳体的组合刚度骨架。若模型把这种高刚度骨架的拱架和侧限作用估计过强，会低估坝顶沉降；组合刚度不足时，软区压缩路径更直接，LOO 预测又可能偏高。"
        }
      ],
      "loo_nrmse": 0.16931250850851604,
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
              "combined_cap_shell_area_fraction": 0.07364588653322478
            },
            "loo_prediction_m": 0.08448614288613403,
            "observed_m": 0.07951806597548103,
            "observed_minus_loo_m": -0.004968076910652999
          },
          "matched_low_residual": {
            "sample_index": 44,
            "design_id": "eval_0044",
            "original_variables": {
              "upstream_slope": 1.3188112255185842,
              "downstream_slope": 1.5114079762250185,
              "crest_width": 6.221145557239652,
              "modulus_zone_height": 15.564859686419368,
              "secondary_top_width": 29.190956093370914,
              "secondary_height": 29.56320457160473,
              "secondary_bottom_width": 33.26097913086414,
              "secondary_right_slope_distance": 4.094621254131198
            },
            "geometric_primitives": {
              "A_tot": 7489.468231786994,
              "A_mod": 392.6549033710563,
              "A_sec": 923.1396684720391,
              "A_drain": 1102.1376076785768,
              "A_main": 5071.536052265321,
              "A_ds_shell": 214.09018419831176,
              "B_bot": 208.25736197599997,
              "B_mod_bot": 47.75904795200003,
              "h_drain": 25.871935742000005,
              "d_ds_perp": 4.094621254131198,
              "d_us_min": 11.509097987000018,
              "X_c_sec": 30.665164438000033,
              "Y_c_sec": 40.33242995400002,
              "W_sec_mean": 31.225967612500014,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "narrow_crest_broad_base_index": 33.475725661754915,
              "soft_area_ratio": 0.12325837294482751,
              "combined_cap_shell_area_fraction": 0.0810131065105804
            },
            "loo_prediction_m": 0.07469771101446823,
            "observed_m": 0.07360604797936852,
            "observed_minus_loo_m": -0.0010916630350997075
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
              "combined_cap_shell_area_fraction": 0.10969660747767776
            },
            "loo_prediction_m": 0.07637088577423387,
            "observed_m": 0.08108230883414239,
            "observed_minus_loo_m": 0.004711423059908515
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
              "combined_cap_shell_area_fraction": 0.09749015503289059
            },
            "loo_prediction_m": 0.07464248772827255,
            "observed_m": 0.07408918199300787,
            "observed_minus_loo_m": -0.0005533057352646725
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
              "combined_cap_shell_area_fraction": 0.10162956256788194
            },
            "loo_prediction_m": 0.08144833974549492,
            "observed_m": 0.07708631363929866,
            "observed_minus_loo_m": -0.004362026106196265
          },
          "matched_low_residual": {
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
              "combined_cap_shell_area_fraction": 0.09946024206847734
            },
            "loo_prediction_m": 0.061430175776745324,
            "observed_m": 0.06226105178888307,
            "observed_minus_loo_m": 0.0008308760121377434
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
              "combined_cap_shell_area_fraction": 0.14036541617054732
            },
            "loo_prediction_m": 0.07977341257497561,
            "observed_m": 0.08383617259362591,
            "observed_minus_loo_m": 0.004062760018650302
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
              "narrow_crest_broad_base_index": 18.07612242309212,
              "soft_area_ratio": 0.1423217374185611,
              "combined_cap_shell_area_fraction": 0.12311484541931035
            },
            "loo_prediction_m": 0.09099510135359834,
            "observed_m": 0.09035788391975826,
            "observed_minus_loo_m": -0.0006372174338400843
          }
        },
        {
          "high_residual": {
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
              "narrow_crest_broad_base_index": 25.457412380429233,
              "soft_area_ratio": 0.09472346962872925,
              "combined_cap_shell_area_fraction": 0.08793725667966058
            },
            "loo_prediction_m": 0.0699775133556377,
            "observed_m": 0.06611362311711586,
            "observed_minus_loo_m": -0.00386389023852185
          },
          "matched_low_residual": {
            "sample_index": 39,
            "design_id": "eval_0039",
            "original_variables": {
              "upstream_slope": 1.3373813543468713,
              "downstream_slope": 1.4245610639452935,
              "crest_width": 7.160981945693493,
              "modulus_zone_height": 14.435939621180296,
              "secondary_top_width": 29.824765343219042,
              "secondary_height": 27.715333439409733,
              "secondary_bottom_width": 26.30811740178615,
              "secondary_right_slope_distance": 5.63574855029583
            },
            "geometric_primitives": {
              "A_tot": 7369.69169958652,
              "A_mod": 344.079878837733,
              "A_sec": 777.8707811233426,
              "A_drain": 993.4735942064738,
              "A_main": 5254.267445418971,
              "A_ds_shell": 266.9912606823638,
              "B_bot": 204.14471639,
              "B_mod_bot": 44.23991639800002,
              "h_drain": 28.848726939000017,
              "d_ds_perp": 5.63574855029583,
              "d_us_min": 4.954312810999994,
              "X_c_sec": 23.585747798,
              "Y_c_sec": 42.995781582999996,
              "W_sec_mean": 28.066441373000032,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "narrow_crest_broad_base_index": 28.507922228845946,
              "soft_area_ratio": 0.1055499758785005,
              "combined_cap_shell_area_fraction": 0.08291678464031016
            },
            "loo_prediction_m": 0.07995754143822065,
            "observed_m": 0.07965216683912185,
            "observed_minus_loo_m": -0.0003053745990987955
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
        },
        {
          "name": "face_asymmetry_soft_elevation_index",
          "expression": "(upstream_slope - downstream_slope) * (Y_c_sec / H_dam)",
          "physical_rationale": "候选假设：五个高残差断面相对各自对照均具有更大的符号化上游坡比减下游坡比，即从较负值变小或转为正值。把这个不对称量按次堆石形心高度缩放，可表征面坡形状在软区所在高程对推力线和偏压路径的影响。正残差的 eval_0050、eval_0025、eval_0051 中，模型可能低估这种偏压导致的整体压缩；负残差的 eval_0031、lhs_009 中，同样的不对称若伴随较强帽盖或骨架约束，预测又可能偏高。"
        }
      ],
      "loo_nrmse": 0.27568490674307244,
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
              "deep_soft_drainage_deficit": 0.10476758365273954,
              "face_asymmetry_soft_elevation_index": -0.0037241751466595167
            },
            "loo_prediction_m": 0.15917534822751672,
            "observed_m": 0.14523452720601138,
            "observed_minus_loo_m": -0.013940821021505345
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
              "deep_soft_drainage_deficit": 0.12187454201233337,
              "face_asymmetry_soft_elevation_index": 3.480005556720041e-05
            },
            "loo_prediction_m": 0.14205490948778784,
            "observed_m": 0.14069650531293582,
            "observed_minus_loo_m": -0.0013584041748520292
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
              "deep_soft_drainage_deficit": 0.16994043175725815,
              "face_asymmetry_soft_elevation_index": 0.04639603369344985
            },
            "loo_prediction_m": 0.1777699623833547,
            "observed_m": 0.1658134202298651,
            "observed_minus_loo_m": -0.011956542153489624
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
              "deep_soft_drainage_deficit": 0.1496438312170615,
              "face_asymmetry_soft_elevation_index": -0.006158995225422185
            },
            "loo_prediction_m": 0.15483197162442922,
            "observed_m": 0.15472294791486899,
            "observed_minus_loo_m": -0.00010902370956022978
          }
        },
        {
          "high_residual": {
            "sample_index": 51,
            "design_id": "eval_0051",
            "original_variables": {
              "upstream_slope": 1.3008948158472777,
              "downstream_slope": 1.3802878566086292,
              "crest_width": 6.430391343310475,
              "modulus_zone_height": 12.455188571475446,
              "secondary_top_width": 19.374227672815323,
              "secondary_height": 32.71979385986924,
              "secondary_bottom_width": 43.14390442799777,
              "secondary_right_slope_distance": 2.147774640470743
            },
            "geometric_primitives": {
              "A_tot": 7097.329771800699,
              "A_mod": 243.61576560880977,
              "A_sec": 1022.7901974150899,
              "A_drain": 1454.3494244889448,
              "A_main": 4376.574384287854,
              "A_ds_shell": 113.0818849913588,
              "B_bot": 197.437904596,
              "B_mod_bot": 36.66165212499999,
              "h_drain": 25.82501756900001,
              "d_ds_perp": 2.147774640470743,
              "d_us_min": 14.037820417000006,
              "X_c_sec": 29.686291771000015,
              "Y_c_sec": 40.11154554900003,
              "W_sec_mean": 31.25906605049998,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_core_eccentricity": 0.15035761158296573,
              "deep_soft_drainage_deficit": 0.1890845802583537,
              "face_asymmetry_soft_elevation_index": -0.043034832037467094
            },
            "loo_prediction_m": 0.15526303044198633,
            "observed_m": 0.16719890679974508,
            "observed_minus_loo_m": 0.011935876357758751
          },
          "matched_low_residual": {
            "sample_index": 33,
            "design_id": "eval_0033",
            "original_variables": {
              "upstream_slope": 1.3706017356365918,
              "downstream_slope": 1.4603382971137762,
              "crest_width": 6.549213381484151,
              "modulus_zone_height": 12.150965901091695,
              "secondary_top_width": 10.416148118674755,
              "secondary_height": 32.25456736981869,
              "secondary_bottom_width": 41.97594594210386,
              "secondary_right_slope_distance": 3.745033367536962
            },
            "geometric_primitives": {
              "A_tot": 7514.709266592214,
              "A_mod": 250.5479703953888,
              "A_sec": 844.9421637691432,
              "A_drain": 1642.3814470602515,
              "A_main": 4776.837685367429,
              "A_ds_shell": 207.46530284159257,
              "B_bot": 208.63877129500003,
              "B_mod_bot": 38.21398451499999,
              "h_drain": 26.594466728999976,
              "d_ds_perp": 3.745033367536962,
              "d_us_min": 21.561961811999993,
              "X_c_sec": 32.086240331,
              "Y_c_sec": 39.48351494100001,
              "W_sec_mean": 26.19604703050001,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_core_eccentricity": 0.15378848395168312,
              "deep_soft_drainage_deficit": 0.15844901140564946,
              "face_asymmetry_soft_elevation_index": -0.04787993061943749
            },
            "loo_prediction_m": 0.14295117865771292,
            "observed_m": 0.14539569678517392,
            "observed_minus_loo_m": 0.002444518127460993
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
              "deep_soft_drainage_deficit": 0.12496067101678475,
              "face_asymmetry_soft_elevation_index": -0.0034164435775396076
            },
            "loo_prediction_m": 0.12023493412941176,
            "observed_m": 0.1320332795531413,
            "observed_minus_loo_m": 0.011798345423729525
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
              "soft_core_eccentricity": 0.10919733731134314,
              "deep_soft_drainage_deficit": 0.1325619957747478,
              "face_asymmetry_soft_elevation_index": -0.017135776344817655
            },
            "loo_prediction_m": 0.14593238383465285,
            "observed_m": 0.1440955601120864,
            "observed_minus_loo_m": -0.0018368237225664386
          }
        },
        {
          "high_residual": {
            "sample_index": 26,
            "design_id": "eval_0026",
            "original_variables": {
              "upstream_slope": 1.304283146560192,
              "downstream_slope": 1.410111927986145,
              "crest_width": 6.392159013077617,
              "modulus_zone_height": 15.078324740752578,
              "secondary_top_width": 18.53951260447502,
              "secondary_height": 28.71927883476019,
              "secondary_bottom_width": 27.295364611782134,
              "secondary_right_slope_distance": 8.086830218322575
            },
            "geometric_primitives": {
              "A_tot": 7185.312514518846,
              "A_mod": 353.93378015664354,
              "A_sec": 658.1723095773623,
              "A_drain": 1239.354886062406,
              "A_main": 4933.851538722435,
              "A_ds_shell": 396.2768821528854,
              "B_bot": 199.85739002000003,
              "B_mod_bot": 44.42751826200001,
              "h_drain": 27.202396424000028,
              "d_ds_perp": 8.086830218322575,
              "d_us_min": 12.272462869999998,
              "X_c_sec": 24.370614429,
              "Y_c_sec": 40.64766025900002,
              "W_sec_mean": 22.917438607999998,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_core_eccentricity": 0.12194002146511168,
              "deep_soft_drainage_deficit": 0.15633138812638409,
              "face_asymmetry_soft_elevation_index": -0.058130977743596
            },
            "loo_prediction_m": 0.15661315717571211,
            "observed_m": 0.14686587884593727,
            "observed_minus_loo_m": -0.009747278329774844
          },
          "matched_low_residual": {
            "sample_index": 27,
            "design_id": "eval_0027",
            "original_variables": {
              "upstream_slope": 1.3056473664939403,
              "downstream_slope": 1.4340507615357638,
              "crest_width": 7.228919478133321,
              "modulus_zone_height": 14.876976828090847,
              "secondary_top_width": 10.501496735960245,
              "secondary_height": 30.423801876604557,
              "secondary_bottom_width": 45.867973673157394,
              "secondary_right_slope_distance": 6.909062474966049
            },
            "geometric_primitives": {
              "A_tot": 7313.767408192272,
              "A_mod": 361.42258919969555,
              "A_sec": 857.4867998220822,
              "A_drain": 1876.491797867288,
              "A_main": 4218.3662213032085,
              "A_ds_shell": 361.80544623935,
              "B_bot": 202.56657644299997,
              "B_mod_bot": 45.156990612000016,
              "h_drain": 25.699221294999973,
              "d_ds_perp": 6.909062474966049,
              "d_us_min": 22.954138701000005,
              "X_c_sec": 27.340317226000025,
              "Y_c_sec": 37.729782288000024,
              "W_sec_mean": 28.184735203999978,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_core_eccentricity": 0.1349695379469144,
              "deep_soft_drainage_deficit": 0.15372722070140307,
              "face_asymmetry_soft_elevation_index": -0.0654680018914603
            },
            "loo_prediction_m": 0.15652241314100676,
            "observed_m": 0.15635562040251327,
            "observed_minus_loo_m": -0.00016679273849348752
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
          "name": "base_spread_crest_throat_index",
          "expression": "(B_bot - crest_width) / (B_bot + crest_width)",
          "physical_rationale": "候选假设：四个匹配对中高残差断面相对对照表现出更明显的宽底窄顶几何，lhs_001 为例外；这在多数负残差断面中与观测位移小于预测相对应。该几何代表坝底展开相对于坝顶受力 throat 的约束尺度，宽底可提供更强的基底抗滑和侧向约束，窄顶则改变上部剪切传递。若模型低估这种整体基底约束，水平位移会被高估；正残差的 lhs_019 提示该效应并非单独决定符号，仍需与软区压缩机制共同解释。"
        },
        {
          "name": "drain_slope_asymmetry_shear",
          "expression": "abs(upstream_slope - downstream_slope) * (A_drain / A_tot)",
          "physical_rationale": "候选假设：正残差 lhs_019、eval_0037、eval_0034 的坡比不对称与底部排水面积占比组合均高于对照，负残差 lhs_017 亦高，负残差 lhs_001 不服从。该无量纲量表征坡比不对称引起的推力偏移与刚性底部排水约束的耦合：排水越强，偏移推力越可能沿底部剪切带转化为下游水平位移；正残差子组中模型可能低估该剪切转化，负残差子组中若主骨架很厚则可能高估。"
        }
      ],
      "loo_nrmse": 0.3337647612526986,
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
            "sample_index": 15,
            "design_id": "lhs_015",
            "original_variables": {
              "upstream_slope": 1.330067888507769,
              "downstream_slope": 1.6538643380389764,
              "crest_width": 10.105264623039737,
              "modulus_zone_height": 8.088284809986758,
              "secondary_top_width": 16.701126741418904,
              "secondary_height": 21.205121572254274,
              "secondary_bottom_width": 55.28494324306234,
              "secondary_right_slope_distance": 6.607914647939461
            },
            "geometric_primitives": {
              "A_tot": 8184.754006807319,
              "A_mod": 156.67531762998547,
              "A_sec": 763.2366827567953,
              "A_drain": 4430.871938240758,
              "A_main": 2833.9700681797785,
              "A_ds_shell": 268.2284303453537,
              "B_bot": 223.51624487799998,
              "B_mod_bot": 31.70107872799997,
              "h_drain": 41.706593618,
              "d_ds_perp": 6.607914647939461,
              "d_us_min": 2.472503617000001,
              "X_c_sec": 11.089266749000018,
              "Y_c_sec": 50.41486541199998,
              "W_sec_mean": 35.9930349925,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_eccentric_area_fraction": 0.004626444343483997,
              "base_spread_crest_throat_index": 0.9134902891037542,
              "drain_slope_asymmetry_shear": 0.17528939791428863
            },
            "loo_prediction_m": 0.06901940242711134,
            "observed_m": 0.04613864520985835,
            "observed_minus_loo_m": -0.022880757217252988
          },
          "matched_low_residual": {
            "sample_index": 8,
            "design_id": "lhs_008",
            "original_variables": {
              "upstream_slope": 1.3023391728933125,
              "downstream_slope": 1.5602052766483099,
              "crest_width": 10.761983864271716,
              "modulus_zone_height": 18.23132658605646,
              "secondary_top_width": 21.964774752174982,
              "secondary_height": 26.670838819588102,
              "secondary_bottom_width": 60.89361800752465,
              "secondary_right_slope_distance": 12.305843009382816
            },
            "geometric_primitives": {
              "A_tot": 7899.8605486849065,
              "A_mod": 620.2211354917332,
              "A_sec": 1104.951419095004,
              "A_drain": 2685.1979186458657,
              "A_main": 3489.4900754523032,
              "A_ds_shell": 604.0308365828677,
              "B_bot": 215.190268621,
              "B_mod_bot": 60.70626413900004,
              "h_drain": 26.097834594000005,
              "d_ds_perp": 12.305843009382816,
              "d_us_min": 16.253014624000002,
              "X_c_sec": 17.237245716000018,
              "Y_c_sec": 37.344821620999994,
              "W_sec_mean": 41.42919638000001,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_eccentric_area_fraction": 0.011203894322494504,
              "base_spread_crest_throat_index": 0.904741079180229,
              "drain_slope_asymmetry_shear": 0.08764984151619051
            },
            "loo_prediction_m": 0.06063147147795552,
            "observed_m": 0.06346261548903362,
            "observed_minus_loo_m": 0.0028311440110780975
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
              "soft_eccentric_area_fraction": 0.016770508210026933,
              "base_spread_crest_throat_index": 0.9231414522776218,
              "drain_slope_asymmetry_shear": 0.04878071252151286
            },
            "loo_prediction_m": 0.10371912948924789,
            "observed_m": 0.12216219021832611,
            "observed_minus_loo_m": 0.018443060729078226
          },
          "matched_low_residual": {
            "sample_index": 27,
            "design_id": "eval_0027",
            "original_variables": {
              "upstream_slope": 1.3056473664939403,
              "downstream_slope": 1.4340507615357638,
              "crest_width": 7.228919478133321,
              "modulus_zone_height": 14.876976828090847,
              "secondary_top_width": 10.501496735960245,
              "secondary_height": 30.423801876604557,
              "secondary_bottom_width": 45.867973673157394,
              "secondary_right_slope_distance": 6.909062474966049
            },
            "geometric_primitives": {
              "A_tot": 7313.767408192272,
              "A_mod": 361.42258919969555,
              "A_sec": 857.4867998220822,
              "A_drain": 1876.491797867288,
              "A_main": 4218.3662213032085,
              "A_ds_shell": 361.80544623935,
              "B_bot": 202.56657644299997,
              "B_mod_bot": 45.156990612000016,
              "h_drain": 25.699221294999973,
              "d_ds_perp": 6.909062474966049,
              "d_us_min": 22.954138701000005,
              "X_c_sec": 27.340317226000025,
              "Y_c_sec": 37.729782288000024,
              "W_sec_mean": 28.184735203999978,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_eccentric_area_fraction": 0.01582421079427936,
              "base_spread_crest_throat_index": 0.9310860374156857,
              "drain_slope_asymmetry_shear": 0.03294443262502514
            },
            "loo_prediction_m": 0.07902276049075438,
            "observed_m": 0.07772136710422413,
            "observed_minus_loo_m": -0.0013013933865302502
          }
        },
        {
          "high_residual": {
            "sample_index": 55,
            "design_id": "eval_0055",
            "original_variables": {
              "upstream_slope": 1.6841140683740377,
              "downstream_slope": 1.3091320749372244,
              "crest_width": 10.039143830537796,
              "modulus_zone_height": 17.781839044764638,
              "secondary_top_width": 28.127651661634445,
              "secondary_height": 20.315555818378925,
              "secondary_bottom_width": 67.16076717711985,
              "secondary_right_slope_distance": 2.1830713348463178
            },
            "geometric_primitives": {
              "A_tot": 8205.497511566498,
              "A_mod": 608.012624574188,
              "A_sec": 967.9185958713975,
              "A_drain": 3374.3284231093803,
              "A_main": 3255.237868011532,
              "A_ds_shell": 70.70336644270763,
              "B_bot": 224.139353926,
              "B_mod_bot": 61.38374015400004,
              "h_drain": 32.90260513700002,
              "d_ds_perp": 2.1830713348463178,
              "d_us_min": 29.892203525000014,
              "X_c_sec": 18.209639337999988,
              "Y_c_sec": 41.673401909000006,
              "W_sec_mean": 47.64420941949997,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_eccentric_area_fraction": 0.009583345463804007,
              "base_spread_crest_throat_index": 0.9142607547087869,
              "drain_slope_asymmetry_shear": 0.15420300802290968
            },
            "loo_prediction_m": 0.09395360324456367,
            "observed_m": 0.11209721788531822,
            "observed_minus_loo_m": 0.01814361464075455
          },
          "matched_low_residual": {
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
              "soft_eccentric_area_fraction": 0.002350994390714315,
              "base_spread_crest_throat_index": 0.9178420125543542,
              "drain_slope_asymmetry_shear": 0.06088190667741409
            },
            "loo_prediction_m": 0.050942595824841896,
            "observed_m": 0.052388155031748075,
            "observed_minus_loo_m": 0.0014455592069061787
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
              "base_spread_crest_throat_index": 0.9126616633274713,
              "drain_slope_asymmetry_shear": 0.05050518280418756
            },
            "loo_prediction_m": 0.09703068914499714,
            "observed_m": 0.08021133624444485,
            "observed_minus_loo_m": -0.016819352900552295
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
              "base_spread_crest_throat_index": 0.8999028611293446,
              "drain_slope_asymmetry_shear": 0.028109840400488077
            },
            "loo_prediction_m": 0.1137043939147526,
            "observed_m": 0.11281189766268739,
            "observed_minus_loo_m": -0.0008924962520652086
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
              "base_spread_crest_throat_index": 0.946473662514429,
              "drain_slope_asymmetry_shear": 0.03593007715387795
            },
            "loo_prediction_m": 0.09501886964948647,
            "observed_m": 0.0807759804523365,
            "observed_minus_loo_m": -0.014242889197149972
          },
          "matched_low_residual": {
            "sample_index": 33,
            "design_id": "eval_0033",
            "original_variables": {
              "upstream_slope": 1.3706017356365918,
              "downstream_slope": 1.4603382971137762,
              "crest_width": 6.549213381484151,
              "modulus_zone_height": 12.150965901091695,
              "secondary_top_width": 10.416148118674755,
              "secondary_height": 32.25456736981869,
              "secondary_bottom_width": 41.97594594210386,
              "secondary_right_slope_distance": 3.745033367536962
            },
            "geometric_primitives": {
              "A_tot": 7514.709266592214,
              "A_mod": 250.5479703953888,
              "A_sec": 844.9421637691432,
              "A_drain": 1642.3814470602515,
              "A_main": 4776.837685367429,
              "A_ds_shell": 207.46530284159257,
              "B_bot": 208.63877129500003,
              "B_mod_bot": 38.21398451499999,
              "h_drain": 26.594466728999976,
              "d_ds_perp": 3.745033367536962,
              "d_us_min": 21.561961811999993,
              "X_c_sec": 32.086240331,
              "Y_c_sec": 39.48351494100001,
              "W_sec_mean": 26.19604703050001,
              "H_dam": 74.0
            },
            "current_feature_values": {
              "soft_eccentric_area_fraction": 0.017291736750294513,
              "base_spread_crest_throat_index": 0.9391303060779131,
              "drain_slope_asymmetry_shear": 0.019612423909507303
            },
            "loo_prediction_m": 0.07925339977032331,
            "observed_m": 0.07615464416718488,
            "observed_minus_loo_m": -0.003098755603138423
          }
        }
      ]
    }
  }
}
```
