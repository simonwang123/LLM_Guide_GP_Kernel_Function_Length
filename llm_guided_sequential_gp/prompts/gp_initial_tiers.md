# 冷启动：为每个变量建议长度尺度先验

## 本轮任务

你将看到由前 20 个有限元样本计算的统计摘要。请针对每个变形响应，分别判断八个几何变量的长度尺度先验中心应选 short、base 还是 long，总共给出 24 个判断。程序将用这些判断训练首次高斯过程模型，并在训练完成后预测第 21 个设计点。

这里的“长度尺度”描述模型在某个变量方向上允许响应变化的快慢。选择 short 表示你认为该方向可能变化较快；选择 long 表示可能变化较平缓；没有足够理由偏离共同基准时选择 base。你选择的是训练前的先验中心，程序仍会根据数据拟合实际长度尺度。

## 本轮提供的数据

下面的 JSON 只包含前 20 个样本的统计摘要及工程背景，不包含逐样本数值：

- design_variables 和 responses 给出八个变量与三个响应的名称及顺序。
- design_bounds_physical 给出每个变量的物理取值范围。
- prior_sign 和 prior_tier 是已有的工程知识。两个二维表均按 design_variables 的顺序排列行、按 responses 的顺序排列列；前者表示预期变化方向，后者表示这一判断的知识强度。它们不直接表示变化快慢。
- statistics 按变量行、响应列给出控制其他七个输入秩后的偏 Spearman 系数，以及单变量二次拟合相对线性拟合的留一预测 R² 增益。null 表示无法计算；这些筛查统计量不能单独确定变量作用。
- initial_observations_sha256 仅供程序核对这 20 个样本，不是分析依据，也不需要写入回答。

{{INITIAL20_CONTEXT_JSON}}

## 判断方法

1. 按 responses 中的顺序，逐一分析每个响应；对该响应，按 design_variables 中的顺序分析八个变量。每个响应—变量关系都必须选择一个类别。
2. 结合变量的工程含义和两组统计线索作判断。没有逐样本观测，不要声称核对过某两个具体样本或局部趋势。
3. 若证据薄弱、统计量为 null 或相互矛盾，选 base。三个响应分别判断，同一变量在不同响应下可以选不同类别。
4. 每条理由用完整句子说明统计或工程依据。不要把相关性写成已经证实的因果关系，也不要预测未经计算的精度提升。

## 输出格式

只输出一个 JSON 对象，不附加说明文字。对象必须恰好包含以下四个字段：

- schema_version：整数 2。
- status：正常完成时为 "OK"；必要输入缺失或相互矛盾、无法完成判断时为 "INVALID_INPUT"。证据不足以偏离 base 属于正常情况。
- assignments：status 为 "OK" 时必须是恰好 24 条记录的数组。每条记录恰好包含 response、variable、prior_category、reason 四个字符串字段。prior_category 只能是 "short"、"base" 或 "long"；reason 必须说明该项判断的依据。
- error：status 为 "OK" 时为 null；status 为 "INVALID_INPUT" 时为具体的错误说明，此时 assignments 必须为空数组。

assignments 的顺序是：先按 responses 的顺序，再按 design_variables 的顺序。24 个响应—变量组合都要出现一次，不能重复或遗漏。输出前核对类别、顺序和记录数量。
