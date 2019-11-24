import pyRserve
import re
from datetime import datetime, timedelta
import numpy

from flask import Flask, jsonify, request
app = Flask(__name__)

epoch = datetime(1970, 1, 1)

@app.route("/anomaly/", methods=['POST'])
def forecast():
    success, payload = get_request_payload()
    if not success: return jsonify(**payload)

    success, conn = get_connection()
    if not success: return jsonify(**conn)

    try:
        response = process_request(payload, conn)
    finally:
        conn.close()

    return jsonify(**response)

# @app.route("/forecast/batches", methods=['POST'])
# def forecast_batches():
#     success, batches = get_request_payload()
#     if not success:
#         return jsonify(**batches)
# 
#     if not isinstance(batches, list):
#         return jsonify(error='Payload needs to be a list!', code=400)
# 
#     success, conn = get_connection()
#     if not success: return jsonify(**conn)
# 
#     try:
#         response = [process_request(one_request, conn) for one_request in batches]
#     finally:
#         conn.close()
# 
#     return jsonify(payload=response, code=200)

def is_valid_date(date):
    return re.match(r'20\d{2}-\d{2}-\d{2}', date)

def get_request_payload():
    try:
        return True, request.get_json()
    except Exception as ex:
        return False, {'error': 'Request payload is invalid json!', 'code': 400, 'details': str(ex)}

def get_connection():
    try:
        return True, pyRserve.connect(host='localhost', port=6311, atomicArray=True)
    except Exception as ex:
        return False, {'error': 'Unable to connect to R server!', 'code': 500, 'details': str(ex)}

def validate_request(req):
    if not isinstance(req, dict):
        return {'error': 'Payload needs to be a dict!', 'code': 400}
    if 'start_date' not in req:
        return {'error': 'start_date is required!', 'code': 400}
    if not is_valid_date(req['start_date']):
        return {'error': 'start_date is in an invalid format, should be YYYY-mm-dd', 'code': 400}
    if 'end_date' not in req:
        return {'error': 'end_date is required!', 'code': 400}
    if not is_valid_date(req['end_date']):
        return {'error': 'end_date is in an invalid format, should be YYYY-mm-dd', 'code': 400}
    if 'series' not in req:
        return {'error': 'series is required!', 'code': 400}
    if not isinstance(req['series'], list):
        return {'error': 'series needs to be a list!', 'code': 400}

    # all good, probably.
    return None

def process_request(req, conn):
    invalid_request = validate_request(req)
    if invalid_request:
        return invalid_request

    result = conn.r.segment_anomaly(
        numpy.array(req['series'], dtype=numpy.float),
        req['start_date'], req['end_date']
    )
    return json.loads(result[0])

#    data = []
#    for i in range(len(result['dates'])):
#        days = result['dates'][i]
#        val = result['y'][i]
#
#        date = (epoch + timedelta(days=days)).strftime("%Y-%m-%d")
#        data.append({'date': date, 'value': val})
#
#    meta = req['meta'] if 'meta' in req else None
#
#    return {'payload': data, 'code': 200, 'meta': meta}


if __name__ == "__main__":
    app.run(host="0.0.0.0", port=6316)

