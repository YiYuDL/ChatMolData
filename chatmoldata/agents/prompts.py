# flake8: noqa
PREFIX = """
You are an expert chemist and your task is to respond to the question or
solve the problem to the best of your ability using the provided tools.
You are able to process image files by using "Image2SMILES" tool.

"""
# When meeting the task about plotting histgram of distribution, use cx_distribution tool.

# FORMAT_INSTRUCTIONS = """
# You can only respond with a single complete
# "Thought, Action, Action Input" format
# OR a single "Final Answer" format.

# Complete format:

# Thought: (reflect on your progress and decide what to do next)
# Action: (the action name, should be one of [{tool_names}])
# Action Input: (the input string to the action)

# OR

# Thouhght: I now know the final answer.
# Final Answer: (the final answer to the original input question)
# """

# FORMAT_INSTRUCTIONS = """
# You can only respond with a single complete
# "Thought1, Action1, Action Input1, ... Thought2, Action2, Action Input2, ..." format (please add number after Thought, Action, Action Input, like Thought1, Action1, Action Input1, starting from 1, add one each time after one complete)
# OR a single "Final Answer" format.

# Complete format:

# Thought1: (reflect on your progress and decide what to do next)
# Action1: (the action name, should be one of [{tool_names}])
# Action Input1: (the input string to the action)
# Thought2: (reflect on your progress and decide what to do next)
# Action2: (the action name, should be one of [{tool_names}])
# Action Input2: (the input string to the action)

# OR

# Thought: (The task is completed)
# Final Answer: (the final answer to the original input question)
# """



FORMAT_INSTRUCTIONS = """
You can only respond with a single complete
"Thought, Action1, Action Input, ... Thought, Action2, Action Input ..." format (please add number after Action, like Action1 starting from 1, add one each time after one complete)
OR a single "Final Answer" format.

Complete format:

Thought: (reflect on your progress and decide what to do next)
Action: (the action name, should be one of [{tool_names}])
Action Input: (the input string to the action)

OR

Thought: The task is completed.
Final Answer: (the final answer to the original input question)
"""

QUESTION_PROMPT = """
Answer the question below using the following tools:

{tool_strings}

Use the tools provided, using the most specific tool available for each action.
Your final answer should contain all information necessary to answer the question and subquestions.

Attention: Importance from front to back, and plan your steps accordingly:
1. when asking for retrieving tasks, "Target2ID" tool should be used and the action input should be the original texts of the question when there is 
no CHEMBL IDs in the content of target.
2. When there is the CHEMBL ID, "SQLtext2CSV" should be used to get the csv file instead of "Target2ID".
3. When Observation is a CSV file and structural images need to be created, the input is the csv file and use the "CSV2Structre_image" tool, return an xlsx file, it should be the last step and stop the task.
4. Before using "CX_prediction" tool, you should use the "CX_preprocess" tool first, then use the file in Observation for prediction task.
5. For substructure search task, you need to use "Subsearch_pre" frist then use "Sub_search".

Question: {input}
"""

SUFFIX = """
Thought: {agent_scratchpad}
"""
FINAL_ANSWER_ACTION = "Final Answer:"


REPHRASE_TEMPLATE = """In this exercise you will assume the role of a scientific assistant. Your task is to answer the provided question as best as you can, based on the provided solution draft.
The solution draft follows the format "Thought, Action, Action Input, Observation", where the 'Thought' statements reflect on your progress and decide what to do next. The rest of the text is information obtained to complement the reasoning sequence, and it is 100% accurate.
Your task is to write an answer to the question based on the solution draft, and the following guidelines:
The text should have an educative and assistant-like tone, be accurate, follow the same reasoning sequence than the solution draft and explain how any conclusion is reached. When input includes "The task execution is completed, please check the error in this step.", please stop the reasoning above.
Question: {question}

Solution draft: {agent_ans}

Answer:
"""
